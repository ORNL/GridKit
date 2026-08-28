#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <utility>

#include <GridKit/Definitions.hpp>
#include <GridKit/Model/VariableMonitorArrowOutput.hpp>

#ifdef GRIDKIT_ENABLE_ARROW

#include <arrow/builder.h>
#include <arrow/io/interfaces.h>
#include <arrow/ipc/writer.h>
#include <arrow/record_batch.h>
#include <arrow/status.h>
#include <arrow/type.h>

namespace GridKit
{
  namespace Model
  {
    namespace
    {
      /// Rows buffered before a record batch is written and flushed
      constexpr int64_t batch_rows = 256;

      /// Convert a failed arrow::Status into an exception
      void throwOnFailure(const arrow::Status& status, const char* what)
      {
        if (!status.ok())
        {
          throw std::runtime_error(std::string(what) + ": " + status.ToString());
        }
      }

      /**
       * @brief Minimal arrow::io::OutputStream over a std::ostream
       *
       * Tracks its own position; the IPC file format footer offsets depend
       * on Tell().
       */
      class OstreamOutputStream : public arrow::io::OutputStream
      {
      public:
        explicit OstreamOutputStream(std::ostream& os)
          : os_(&os)
        {
        }

        arrow::Status Write(const void* data, int64_t nbytes) override
        {
          os_->write(static_cast<const char*>(data), static_cast<std::streamsize>(nbytes));
          if (os_->fail())
          {
            return arrow::Status::IOError("Failed writing to Arrow monitor output stream");
          }
          position_ += nbytes;
          return arrow::Status::OK();
        }

        arrow::Status Flush() override
        {
          os_->flush();
          if (os_->fail())
          {
            return arrow::Status::IOError("Failed flushing Arrow monitor output stream");
          }
          return arrow::Status::OK();
        }

        arrow::Status Close() override
        {
          if (closed_)
          {
            return arrow::Status::OK();
          }

          auto status = Flush();
          if (!status.ok())
          {
            return status;
          }

          closed_ = true;
          return arrow::Status::OK();
        }

        arrow::Result<int64_t> Tell() const override
        {
          return position_;
        }

        bool closed() const override
        {
          return closed_;
        }

      private:
        /// Output stream written to (owned by VariableMonitorArrowOutputImpl)
        std::ostream* os_;
        /// Bytes written so far
        int64_t       position_{0};
        /// Has Close() been called?
        bool          closed_{false};
      };
    } // namespace

    /**
     * @brief Arrow-dependent state and logic behind VariableMonitorArrowOutput
     */
    class VariableMonitorArrowOutputImpl
    {
    public:
      VariableMonitorArrowOutputImpl(const std::string& file_name, bool stream_format)
        : file_(file_name, std::ios::binary | std::ios::trunc),
          stream_format_(stream_format)
      {
        if (!file_.is_open())
        {
          throw std::runtime_error("Failed to open Arrow monitor output file: " + file_name);
        }
        out_ = std::make_shared<OstreamOutputStream>(file_);
      }

      /**
       * @brief Build the all-float64 schema and create the IPC writer
       */
      void beginTable(const std::vector<std::string>& names)
      {
        arrow::FieldVector fields;
        fields.reserve(names.size());
        for (const auto& name : names)
        {
          fields.push_back(arrow::field(name, arrow::float64()));
        }
        schema_ = arrow::schema(std::move(fields));

        auto writer = stream_format_ ? arrow::ipc::MakeStreamWriter(out_, schema_)
                                     : arrow::ipc::MakeFileWriter(out_, schema_);
        throwOnFailure(writer.status(), "Failed to create Arrow IPC writer");
        writer_ = std::move(writer).ValueOrDie();

        columns_.resize(names.size());
        for (auto& column : columns_)
        {
          column.reserve(static_cast<size_t>(batch_rows));
        }
      }

      /**
       * @brief Buffer one row; write a record batch every `batch_rows` rows
       */
      void appendRow(const std::vector<double>& values)
      {
        if (!writer_)
        {
          throw std::runtime_error("Arrow monitor output used before the header was written");
        }
        if (values.size() != columns_.size())
        {
          throw std::runtime_error("Arrow monitor row width does not match the schema");
        }

        for (size_t i = 0; i < values.size(); ++i)
        {
          columns_[i].push_back(values[i]);
        }

        if (static_cast<int64_t>(columns_.front().size()) >= batch_rows)
        {
          flushBatch();
        }
      }

      /**
       * @brief Flush buffered rows and finalize the IPC output
       *
       * The writer Close() emits the file footer (IPC file format) or the
       * end-of-stream marker (IPC stream format) and must precede closing
       * the stream.
       */
      void finish()
      {
        if (writer_)
        {
          flushBatch();
          throwOnFailure(writer_->Close(), "Failed to finalize Arrow IPC output");
          writer_.reset();
        }
        if (out_)
        {
          throwOnFailure(out_->Close(), "Failed to close Arrow monitor output stream");
          out_.reset();
        }
        file_.close();
      }

    private:
      /**
       * @brief Write buffered rows as one record batch and flush the stream
       *
       * The flush makes each batch immediately visible to live consumers
       * reading from a FIFO.
       */
      void flushBatch()
      {
        const auto rows = static_cast<int64_t>(columns_.front().size());
        if (rows == 0)
        {
          return;
        }

        arrow::ArrayVector arrays;
        arrays.reserve(columns_.size());
        for (auto& column : columns_)
        {
          arrow::DoubleBuilder builder;
          throwOnFailure(builder.AppendValues(column), "Failed to append Arrow column values");
          std::shared_ptr<arrow::Array> array;
          throwOnFailure(builder.Finish(&array), "Failed to build Arrow column array");
          arrays.push_back(std::move(array));
          column.clear();
        }

        auto batch = arrow::RecordBatch::Make(schema_, rows, std::move(arrays));
        throwOnFailure(writer_->WriteRecordBatch(*batch), "Failed to write Arrow record batch");
        throwOnFailure(out_->Flush(), "Failed to flush Arrow record batch");
      }

      /// Output file stream (opened in binary mode)
      std::ofstream file_;

      /// Write the IPC stream format instead of the IPC file format
      bool stream_format_{false};

      /// Arrow view of the output file stream
      std::shared_ptr<OstreamOutputStream> out_;

      /// Schema shared by all record batches (t first, then one float64
      /// field per monitored variable)
      std::shared_ptr<arrow::Schema> schema_;

      /// IPC writer for the selected format
      std::shared_ptr<arrow::ipc::RecordBatchWriter> writer_;

      /// Row buffer, one vector per column
      std::vector<std::vector<double>> columns_;
    };

    VariableMonitorArrowOutput::VariableMonitorArrowOutput(const std::string& file_name, bool stream_format)
      : file_name_(file_name),
        stream_format_(stream_format)
    {
    }

    VariableMonitorArrowOutput::VariableMonitorArrowOutput(VariableMonitorArrowOutput&&) noexcept = default;

    VariableMonitorArrowOutput::~VariableMonitorArrowOutput() = default;

    void VariableMonitorArrowOutput::start()
    {
      if (!impl_)
      {
        impl_ = std::make_unique<VariableMonitorArrowOutputImpl>(file_name_, stream_format_);
      }
    }

    void VariableMonitorArrowOutput::beginTable(const std::vector<std::string>& names) const
    {
      if (!impl_)
      {
        throw std::runtime_error("Arrow monitor output used before start()");
      }
      impl_->beginTable(names);
    }

    void VariableMonitorArrowOutput::appendRow(const std::vector<double>& values) const
    {
      if (!impl_)
      {
        throw std::runtime_error("Arrow monitor output used before start()");
      }
      impl_->appendRow(values);
    }

    void VariableMonitorArrowOutput::stop()
    {
      if (impl_)
      {
        impl_->finish();
        impl_.reset();
      }
    }
  } // namespace Model
} // namespace GridKit

#else // GRIDKIT_ENABLE_ARROW

namespace GridKit
{
  namespace Model
  {
    /**
     * @brief Placeholder giving unique_ptr a complete type; never
     * instantiated because construction of VariableMonitorArrowOutput throws
     */
    class VariableMonitorArrowOutputImpl
    {
    };

    VariableMonitorArrowOutput::VariableMonitorArrowOutput(const std::string&, bool)
    {
      throw std::runtime_error(
          "GridKit was built without Apache Arrow support; "
          "rebuild with GridKit_ENABLE_ARROW=ON to use Arrow monitor output");
    }

    VariableMonitorArrowOutput::VariableMonitorArrowOutput(VariableMonitorArrowOutput&&) noexcept = default;

    VariableMonitorArrowOutput::~VariableMonitorArrowOutput() = default;

    void VariableMonitorArrowOutput::start()
    {
    }

    void VariableMonitorArrowOutput::beginTable(const std::vector<std::string>&) const
    {
    }

    void VariableMonitorArrowOutput::appendRow(const std::vector<double>&) const
    {
    }

    void VariableMonitorArrowOutput::stop()
    {
    }
  } // namespace Model
} // namespace GridKit

#endif // GRIDKIT_ENABLE_ARROW
