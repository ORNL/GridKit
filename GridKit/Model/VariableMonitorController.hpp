#pragma once

/**
 * @file
 */

#include <ostream>
#include <string>

namespace GridKit
{
  namespace Model
  {
    /**
     * @brief Monitor associated with the SystemModel; controls output from
     * component and bus monitors.
     *
     * High-level print functions (without parameters) manage printing for all
     * monitors for multiple output sinks.
     */
    template <typename scalar_type>
    class VariableMonitorController : public VariableMonitorBase
    {
    public:
      /// Underlying scalar value type
      using ScalarT  = scalar_type;
      /// Underlying real value type
      using RealT    = typename GridKit::ScalarTraits<ScalarT>::RealT;
      ///@{
      /// @brief Alias
      using Format   = VariableMonitorFormat;
      using Layout   = VariableMonitorLayout;
      using SinkSpec = VariableMonitorBase::SinkSpec;
      ///@}

      /// Default to empty monitor
      VariableMonitorController() = default;

      /**
       * @brief Constructor expects a time variable to monitor
       */
      explicit VariableMonitorController(const RealT& time_var)
        : time_(&time_var)
      {
      }

      virtual ~VariableMonitorController()
      {
      }

      /**
       * @brief Add a monitor to the output
       *
       * Each component and bus could have a monitor for their respective
       * values.
       *
       * @param monitor Monitor to add (raw pointer indicates ownership is
       * elsewhere)
       */
      void addMonitor(const VariableMonitorBase* monitor)
      {
        if (monitor && !monitor->empty())
        {
          monitor->appendEntries(entries_);
        }
      }

      /**
       * @brief Add output sink based on spec
       *
       * If `spec.file_name` is empty, `std::cout` is used.
       *
       * @param spec Specifies details for the sink.
       * @param os Output stream used when `spec.file_name` is empty.
       */
      void addSink(const SinkSpec& spec, std::ostream& os = std::cout)
      {
        if (spec.file_name.empty())
        {
          sinks_.push_back(make_sink(spec, os));
        }
        else
        {
          sinks_.push_back(make_sink(spec, spec.file_name));
        }
      }

      bool empty() const override
      {
        return (!time_) || entries_.empty();
      }

      /**
       * @brief Print header if we're monitoring
       */
      void start()
      {
        if (!empty())
        {
          startSinks();
          printHeader();
        }
      }

      /**
       * @brief Print footer if we're monitoring
       */
      void stop()
      {
        if (!empty())
        {
          printFooter();
          stopSinks();
        }
      }

      /**
       * @brief Organize header output for this and all submonitors
       */
      template <typename SinkT>
      void printFullHeader(const SinkT& sink) const
      {
        buffer_.clear();
        appendHeader(buffer_, sink);
        buffer_ += '\n';
        sink.os->write(buffer_.data(), static_cast<std::streamsize>(buffer_.size()));
      }

      /**
       * @brief Print header for all sinks
       */
      void printHeader() const
      {
        for (auto&& sink : sinks_)
        {
          std::visit([this](auto&& sink)
                     { this->printFullHeader(sink); },
                     sink);
        }
      }

      /**
       * @brief Organize variable output for this and all submonitors
       */
      template <typename SinkT>
      void printFull(const SinkT& sink) const
      {
        buffer_.clear();
        if constexpr (SinkT::layout == Layout::LONG)
        {
          const auto  sample = std::to_string(sample_);
          std::string t;
          VariableMonitorDetail::formatReal(t, *time_);
          for (auto&& entry : entries_)
          {
            appendRecord(buffer_, sink, sample, t, entry);
          }
        }
        else
        {
          append(buffer_, sink);
          buffer_ += '\n';
        }
        sink.os->write(buffer_.data(), static_cast<std::streamsize>(buffer_.size()));
      }

      /**
       * @brief Print variables to each sink
       */
      void print() const
      {
        for (auto&& sink : sinks_)
        {
          std::visit([this](auto&& sink)
                     { this->printFull(sink); },
                     sink);
        }
        ++sample_;
      }

      /**
       * @brief Organize footer output for this and all submonitors
       */
      template <typename SinkT>
      void printFullFooter(const SinkT& sink) const
      {
        buffer_.clear();
        appendFooter(buffer_, sink);
        sink.os->write(buffer_.data(), static_cast<std::streamsize>(buffer_.size()));
      }

      /**
       * @brief Print footer for all sinks
       */
      void printFooter() const
      {
        for (auto&& sink : sinks_)
        {
          std::visit([this](auto&& sink)
                     { this->printFullFooter(sink); },
                     sink);
        }
      }

    protected:
      /**
       * @brief Append column names (CSV) or the opening bracket (JSON)
       */
      template <typename SinkT>
      void appendHeader(std::string& out, const SinkT& sink) const
      {
        if constexpr (SinkT::format == Format::CSV && SinkT::layout == Layout::WIDE)
        {
          out += "t";
          for (auto&& entry : entries_)
          {
            out += sink.delim;
            out += entry.device_class;
            out += '_';
            out += entry.id;
            out += '_';
            out += entry.variable;
          }
        }
        else if constexpr (SinkT::format == Format::CSV)
        {
          for (auto* column : {"sample", "t", "class", "id", "variable"})
          {
            out += column;
            out += sink.delim;
          }
          out += "value";
        }
        else if constexpr (SinkT::format == Format::JSON)
        {
          out += "[\n";
        }
      }

      /**
       * @brief Append one wide layout row for the current sample
       */
      template <typename SinkT>
      void append(std::string& out, const SinkT& sink) const
      {
        if constexpr (SinkT::format == Format::CSV)
        {
          VariableMonitorDetail::formatReal(out, *time_);
          for (auto&& entry : entries_)
          {
            out += sink.delim;
            VariableMonitorDetail::formatReal(out, entry.value());
          }
        }
        else if constexpr (SinkT::format == Format::JSON)
        {
          if (sample_ > 0)
          {
            out += "  ,\n";
          }
          out += "  {\n    \"t\": ";
          VariableMonitorDetail::formatReal(out, *time_);

          const Entry* previous = nullptr;
          for (auto&& entry : entries_)
          {
            if (previous && previous->device_class == entry.device_class && previous->id == entry.id)
            {
              out += ",\n      \"";
            }
            else
            {
              if (previous)
              {
                out += "\n    }";
              }
              out += ",\n    \"";
              out += entry.device_class;
              out += '_';
              out += entry.id;
              out += "\": {\n      \"";
            }
            out += entry.variable;
            out += "\": ";
            VariableMonitorDetail::formatReal(out, entry.value());

            previous = &entry;
          }
          if (previous)
          {
            out += "\n    }";
          }
          out += "\n  }";
        }
        else
        {
          out += "  - t: ";
          VariableMonitorDetail::formatReal(out, *time_);
          out += '\n';

          const Entry* previous = nullptr;
          for (auto&& entry : entries_)
          {
            if (!previous || previous->device_class != entry.device_class || previous->id != entry.id)
            {
              out += "    ";
              out += entry.device_class;
              out += '_';
              out += entry.id;
              out += ":\n";
            }
            out += "      ";
            out += entry.variable;
            out += ": ";
            VariableMonitorDetail::formatReal(out, entry.value());
            out += '\n';

            previous = &entry;
          }
        }
      }

      /**
       * @brief Append one record of the long layout
       */
      template <typename SinkT>
      void appendRecord(std::string& out, const SinkT& sink, const std::string& sample, const std::string& t, const Entry& entry) const
      {
        if constexpr (SinkT::format == Format::CSV)
        {
          out += sample;
          out += sink.delim;
          out += t;
          out += sink.delim;
          out += entry.device_class;
          out += sink.delim;
          out += entry.id;
          out += sink.delim;
          out += entry.variable;
          out += sink.delim;
          VariableMonitorDetail::formatReal(out, entry.value());
          out += '\n';
        }
        else if constexpr (SinkT::format == Format::JSON)
        {
          // Comma before every record but the first
          if (sample_ > 0 || &entry != &entries_.front())
          {
            out += ",\n";
          }
          out += "  {\"sample\": ";
          out += sample;
          out += ", \"t\": ";
          out += t;
          out += ", \"class\": \"";
          out += entry.device_class;
          out += "\", \"id\": \"";
          out += entry.id;
          out += "\", \"variable\": \"";
          out += entry.variable;
          out += "\", \"value\": ";
          VariableMonitorDetail::formatReal(out, entry.value());
          out += '}';
        }
        else
        {
          out += "- {sample: ";
          out += sample;
          out += ", t: ";
          out += t;
          out += ", class: ";
          out += entry.device_class;
          out += ", id: \"";
          out += entry.id;
          out += "\", variable: ";
          out += entry.variable;
          out += ", value: ";
          VariableMonitorDetail::formatReal(out, entry.value());
          out += "}\n";
        }
      }

      /**
       * @brief Append the closing bracket (JSON)
       */
      template <typename SinkT>
      void appendFooter(std::string& out, const SinkT&) const
      {
        if constexpr (SinkT::format == Format::JSON)
        {
          out += "\n]\n";
        }
      }

      void appendEntries(std::vector<Entry>& entries) const override
      {
        entries.insert(entries.end(), entries_.begin(), entries_.end());
      }

    private:
      /// Time variable; printed first
      const RealT* time_{nullptr};

      /**
       * @brief Define sink for a specific output format and layout
       */
      template <Format format_type, Layout layout_type>
      struct Sink
      {
        /// Output format
        static constexpr Format format = format_type;
        /// Output layout
        static constexpr Layout layout = layout_type;

        Sink() = delete;

        /**
         * @brief Version for an output stream that already exists
         */
        Sink(std::ostream& out, const std::string& delimiter)
          : os(&out), delim(delimiter)
        {
        }

        /**
         * @brief Version for opening an output stream for the given file
         */
        Sink(const std::string& fileName, const std::string& delimiter)
          : file_name(fileName),
            delim(delimiter)
        {
        }

        /**
         * @brief Have to move because of unique_ptr
         */
        Sink(Sink&&) = default;

        void start()
        {
          if (file_name.empty())
          {
            return;
          }
          if (file_stream)
          {
            return;
          }
          file_stream = std::make_unique<std::ofstream>(file_name);
          os          = file_stream.get();
        }

        void stop()
        {
          if (file_name.empty())
          {
            return;
          }
          file_stream.reset();
          os = nullptr;
        }

        /// Output file name
        std::string                    file_name;
        /// Output file stream (if we opened one)
        std::unique_ptr<std::ofstream> file_stream;
        /// Output stream for printing
        std::ostream*                  os{nullptr};
        /// Delimiter for CSV output
        std::string                    delim;
      };

      /// Variant type for all possible sink types
      using SinkVariant = std::variant<Sink<Format::CSV, Layout::WIDE>,
                                       Sink<Format::JSON, Layout::WIDE>,
                                       Sink<Format::YAML, Layout::WIDE>,
                                       Sink<Format::CSV, Layout::LONG>,
                                       Sink<Format::JSON, Layout::LONG>,
                                       Sink<Format::YAML, Layout::LONG>>;

      /**
       * @brief Construct the sink type whose format and layout match the spec
       */
      template <std::size_t index = 0, typename T>
      static SinkVariant make_sink(const SinkSpec& spec, T&& arg)
      {
        using SinkT = std::variant_alternative_t<index, SinkVariant>;
        if (spec.format == SinkT::format && spec.layout == SinkT::layout)
        {
          return SinkT(std::forward<T>(arg), spec.delim);
        }
        if constexpr (index + 1 < std::variant_size_v<SinkVariant>)
        {
          return make_sink<index + 1>(spec, std::forward<T>(arg));
        }
        throw std::invalid_argument("Invalid monitor output format or layout");
      }

      /**
       * @brief Get sinks ready for output
       */
      void startSinks()
      {
        for (auto&& sink : sinks_)
        {
          std::visit([this](auto&& sink)
                     { (void) this; sink.start(); },
                     sink);
        }
      }

      /**
       * @brief Finalize sinks and close owned files
       */
      void stopSinks()
      {
        for (auto&& sink : sinks_)
        {
          std::visit([this](auto&& sink)
                     { (void) this; sink.stop(); },
                     sink);
        }
      }

      /// Collection of output sinks
      std::vector<SinkVariant> sinks_;
      /// Monitored variables of all submonitors
      std::vector<Entry>       entries_;

      /// Sample index; event times repeat t
      mutable std::size_t sample_{0};

      /// Reused output buffer
      mutable std::string buffer_;
    };
  } // namespace Model
} // namespace GridKit
