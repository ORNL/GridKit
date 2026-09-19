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
      VariableMonitorArrowOutputImpl(const std::string& file_name, bool stream_format, VariableMonitorLayout layout, int64_t batch_rows)
        : file_(file_name, std::ios::binary | std::ios::trunc),
          stream_format_(stream_format),
          layout_(layout),
          batch_rows_(batch_rows)
      {
        if (!file_.is_open())
        {
          throw std::runtime_error("Failed to open Arrow monitor output file: " + file_name);
        }
        out_ = std::make_shared<OstreamOutputStream>(file_);
      }

      /// Build the selected table schema and create the IPC writer
      void beginTable(std::span<const VariableMonitorBase::Entry> entries)
      {
        arrow::FieldVector fields;
        if (layout_ == VariableMonitorLayout::LONG)
        {
          fields = {arrow::field("sample", arrow::uint64()),
                    arrow::field("t", arrow::float64()),
                    arrow::field("class", arrow::utf8()),
                    arrow::field("id", arrow::utf8()),
                    arrow::field("variable", arrow::utf8()),
                    arrow::field("value", arrow::float64())};
        }
        else
        {
          fields.push_back(arrow::field("t", arrow::float64()));
          for (const auto& entry : entries)
          {
            fields.push_back(arrow::field(std::string(entry.device_class) + "_"
                                              + std::string(entry.id) + "_" + std::string(entry.variable),
                                          arrow::float64()));
          }
        }
        schema_     = arrow::schema(std::move(fields));
        auto writer = stream_format_ ? arrow::ipc::MakeStreamWriter(out_, schema_)
                                     : arrow::ipc::MakeFileWriter(out_, schema_);
        throwOnFailure(writer.status(), "Failed to create Arrow IPC writer");
        writer_ = std::move(writer).ValueOrDie();
        for (const auto& field : schema_->fields())
        {
          auto builder = arrow::MakeBuilder(field->type());
          throwOnFailure(builder.status(), "Failed to create Arrow column builder");
          columns_.push_back(std::move(builder).ValueOrDie());
        }
      }

      /// Buffer one sample; long layout can span multiple record batches
      void append(std::size_t sample, double t, std::span<const VariableMonitorBase::Entry> entries)
      {
        if (!writer_)
        {
          throw std::runtime_error("Arrow monitor output used before the header was written");
        }
        if (layout_ == VariableMonitorLayout::LONG)
        {
          for (const auto& entry : entries)
          {
            appendValue<arrow::UInt64Builder>(0, static_cast<uint64_t>(sample));
            appendValue<arrow::DoubleBuilder>(1, t);
            appendValue<arrow::StringBuilder>(2, entry.device_class);
            appendValue<arrow::StringBuilder>(3, entry.id);
            appendValue<arrow::StringBuilder>(4, entry.variable);
            appendValue<arrow::DoubleBuilder>(5, entry.value());
            finishRow();
          }
        }
        else
        {
          if (entries.size() + 1 != columns_.size())
          {
            throw std::runtime_error("Arrow monitor row width does not match the schema");
          }
          appendValue<arrow::DoubleBuilder>(0, t);
          for (std::size_t i = 0; i < entries.size(); ++i)
          {
            appendValue<arrow::DoubleBuilder>(i + 1, entries[i].value());
          }
          finishRow();
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
        if (file_.fail())
        {
          throw std::runtime_error("Failed to close Arrow monitor output file");
        }
      }

    private:
      template <typename BuilderT, typename ValueT>
      void appendValue(std::size_t column, ValueT value)
      {
        throwOnFailure(static_cast<BuilderT&>(*columns_[column]).Append(value),
                       "Failed to append Arrow column value");
      }

      void finishRow()
      {
        if (++rows_ >= batch_rows_)
        {
          flushBatch();
        }
      }

      /**
       * @brief Write buffered rows as one record batch and flush the stream
       *
       * The flush makes each batch immediately visible to live consumers
       * reading from a FIFO.
       */
      void flushBatch()
      {
        if (rows_ == 0)
        {
          return;
        }

        arrow::ArrayVector arrays;
        arrays.reserve(columns_.size());
        for (auto& column : columns_)
        {
          std::shared_ptr<arrow::Array> array;
          throwOnFailure(column->Finish(&array), "Failed to build Arrow column array");
          arrays.push_back(std::move(array));
        }

        auto batch = arrow::RecordBatch::Make(schema_, rows_, std::move(arrays));
        throwOnFailure(writer_->WriteRecordBatch(*batch), "Failed to write Arrow record batch");
        throwOnFailure(out_->Flush(), "Failed to flush Arrow record batch");
        rows_ = 0;
      }

      /// Output file stream (opened in binary mode)
      std::ofstream file_;

      /// Write the IPC stream format instead of the IPC file format
      bool stream_format_{false};

      /// Table layout
      VariableMonitorLayout layout_;

      /// Rows buffered before a record batch is written and flushed
      int64_t batch_rows_{256};

      /// Arrow view of the output file stream
      std::shared_ptr<OstreamOutputStream> out_;

      /// Schema shared by all record batches
      std::shared_ptr<arrow::Schema> schema_;

      /// IPC writer for the selected format
      std::shared_ptr<arrow::ipc::RecordBatchWriter> writer_;

      /// Typed column buffers for the selected schema
      std::vector<std::unique_ptr<arrow::ArrayBuilder>> columns_;
      int64_t                                           rows_{0};
    };

    VariableMonitorArrowOutput::VariableMonitorArrowOutput(const std::string&    file_name,
                                                           bool                  stream_format,
                                                           VariableMonitorLayout layout,
                                                           int64_t               batch_rows)
      : file_name_(file_name),
        stream_format_(stream_format),
        layout_(layout),
        batch_rows_(batch_rows)
    {
      if (file_name.empty() || batch_rows < 1
          || (layout != VariableMonitorLayout::WIDE && layout != VariableMonitorLayout::LONG))
      {
        throw std::invalid_argument("Arrow monitor requires a file name, positive batch_rows and a valid layout");
      }
    }

    VariableMonitorArrowOutput::VariableMonitorArrowOutput(VariableMonitorArrowOutput&&) noexcept = default;

    VariableMonitorArrowOutput::~VariableMonitorArrowOutput() = default;

    void VariableMonitorArrowOutput::start()
    {
      if (!impl_)
      {
        impl_ = std::make_unique<VariableMonitorArrowOutputImpl>(file_name_, stream_format_, layout_, batch_rows_);
      }
    }

    void VariableMonitorArrowOutput::beginTable(std::span<const VariableMonitorBase::Entry> entries) const
    {
      if (!impl_)
      {
        throw std::runtime_error("Arrow monitor output used before start()");
      }
      impl_->beginTable(entries);
    }

    void VariableMonitorArrowOutput::append(std::size_t sample, double t, std::span<const VariableMonitorBase::Entry> entries) const
    {
      if (!impl_)
      {
        throw std::runtime_error("Arrow monitor output used before start()");
      }
      impl_->append(sample, t, entries);
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

    VariableMonitorArrowOutput::VariableMonitorArrowOutput(const std::string&, bool, VariableMonitorLayout, int64_t)
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

    void VariableMonitorArrowOutput::beginTable(std::span<const VariableMonitorBase::Entry>) const
    {
    }

    void VariableMonitorArrowOutput::append(std::size_t, double, std::span<const VariableMonitorBase::Entry>) const
    {
    }

    void VariableMonitorArrowOutput::stop()
    {
    }
  } // namespace Model
} // namespace GridKit

#endif // GRIDKIT_ENABLE_ARROW
