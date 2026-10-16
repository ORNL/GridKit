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
      using Csv      = VariableMonitorBase::Csv;
      using Json     = VariableMonitorBase::Json;
      using Yaml     = VariableMonitorBase::Yaml;
      using Wide     = VariableMonitorBase::Wide;
      using Long     = VariableMonitorBase::Long;
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
      template <typename FormatT>
      void printFullHeader(std::ostream& os, const FormatT& fmt, Wide) const
      {
        buffer_.clear();
        appendHeader(buffer_, fmt);
        buffer_ += '\n';
        os.write(buffer_.data(), static_cast<std::streamsize>(buffer_.size()));
      }

      /**
       * @brief Organize long layout header output; only CSV adds column names
       */
      template <typename FormatT>
      void printFullHeader(std::ostream& os, const FormatT& fmt, Long) const
      {
        printFullHeader(os, fmt, Wide{});
      }

      void printFullHeader(std::ostream& os, const Csv& csv, Long) const
      {
        buffer_.clear();
        for (auto* column : {"sample", "t", "class", "id", "variable"})
        {
          buffer_ += column;
          buffer_ += csv.delim;
        }
        buffer_ += "value\n";
        os.write(buffer_.data(), static_cast<std::streamsize>(buffer_.size()));
      }

      /**
       * @brief Print header for all sinks
       */
      void printHeader() const
      {
        for (auto&& sink : sinks_)
        {
          std::visit([this](auto&& sink)
                     { this->printFullHeader(*sink.os, sink.format, sink.layout); },
                     sink);
        }
      }

      /**
       * @brief Organize variable output for this and all submonitors
       */
      template <typename FormatT>
      void printFull(std::ostream& os, const FormatT& fmt, Wide) const
      {
        buffer_.clear();
        append(buffer_, fmt);
        buffer_ += '\n';
        os.write(buffer_.data(), static_cast<std::streamsize>(buffer_.size()));
      }

      /**
       * @brief Organize long layout variable output, one record per variable
       */
      template <typename FormatT>
      void printFull(std::ostream& os, const FormatT& fmt, Long) const
      {
        const auto  sample = std::to_string(sample_);
        std::string t;
        VariableMonitorDetail::formatReal(t, *time_);
        buffer_.clear();
        for (auto&& entry : entries_)
        {
          appendRecord(buffer_, fmt, sample, t, entry);
        }
        os.write(buffer_.data(), static_cast<std::streamsize>(buffer_.size()));
      }

      /**
       * @brief Print variables to each sink
       */
      void print() const
      {
        for (auto&& sink : sinks_)
        {
          std::visit([this](auto&& sink)
                     {
              this->printFull(*sink.os, sink.format, sink.layout);
              using T = std::remove_cvref_t<decltype(sink)>;
              if constexpr (std::is_same_v<T, Sink<Json, Wide>>)
              {
                sink.format.after_first = true;
              } },
                     sink);
        }
        ++sample_;
      }

      /**
       * @brief Organize footer output for this and all submonitors
       */
      template <typename FormatT>
      void printFullFooter(std::ostream& os, const FormatT& fmt, Wide) const
      {
        buffer_.clear();
        appendFooter(buffer_, fmt);
        os.write(buffer_.data(), static_cast<std::streamsize>(buffer_.size()));
      }

      /**
       * @brief Organize long layout footer output, the same as wide
       */
      template <typename FormatT>
      void printFullFooter(std::ostream& os, const FormatT& fmt, Long) const
      {
        printFullFooter(os, fmt, Wide{});
      }

      /**
       * @brief Print footer for all sinks
       */
      void printFooter() const
      {
        for (auto&& sink : sinks_)
        {
          std::visit([this](auto&& sink)
                     { this->printFullFooter(*sink.os, sink.format, sink.layout); },
                     sink);
        }
      }

    protected:
      void appendHeader(std::string& out, Csv csv) const
      {
        out += "t";
        for (auto&& entry : entries_)
        {
          out += csv.delim;
          out += entry.device_class;
          out += '_';
          out += entry.id;
          out += '_';
          out += entry.variable;
        }
      }

      void appendHeader(std::string& out, Json) const
      {
        out += "[\n";
      }

      void appendHeader(std::string&, Yaml) const
      {
      }

      void append(std::string& out, Csv csv) const
      {
        VariableMonitorDetail::formatReal(out, *time_);
        for (auto&& entry : entries_)
        {
          out += csv.delim;
          VariableMonitorDetail::formatReal(out, entry.value());
        }
      }

      void append(std::string& out, Json json) const
      {
        if (json.after_first)
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

      void append(std::string& out, Yaml) const
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

      void appendEntries(std::vector<Entry>& entries) const override
      {
        entries.insert(entries.end(), entries_.begin(), entries_.end());
      }

      ///@{
      /**
       * @brief Append one record of the long layout
       */
      void appendRecord(std::string& out, const Csv& csv, const std::string& sample, const std::string& t, const Entry& entry) const
      {
        out += sample;
        out += csv.delim;
        out += t;
        out += csv.delim;
        out += entry.device_class;
        out += csv.delim;
        out += entry.id;
        out += csv.delim;
        out += entry.variable;
        out += csv.delim;
        VariableMonitorDetail::formatReal(out, entry.value());
        out += '\n';
      }

      void appendRecord(std::string& out, const Json& json, const std::string& sample, const std::string& t, const Entry& entry) const
      {
        if (json.after_first)
        {
          out += ",\n";
        }
        json.after_first = true;

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

      void appendRecord(std::string& out, const Yaml&, const std::string& sample, const std::string& t, const Entry& entry) const
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

      ///@}

      void appendFooter(std::string&, Csv) const
      {
      }

      void appendFooter(std::string& out, Json) const
      {
        out += "\n]\n";
      }

      void appendFooter(std::string&, Yaml) const
      {
      }

    private:
      /// Time variable; printed first
      const RealT* time_{nullptr};

      /**
       * @brief Define sink for a specific output format and layout
       */
      template <typename format_type, typename layout_type>
      struct Sink
      {
        using FormatT = format_type;
        using LayoutT = layout_type;

        Sink() = delete;

        /**
         * @brief Version for an output stream that already exists
         */
        Sink(std::ostream& out, FormatT fmt, LayoutT layout)
          : os(&out), format(fmt), layout(layout)
        {
        }

        /**
         * @brief Version for opening an output stream for the given file
         */
        Sink(const std::string& fileName, FormatT fmt, LayoutT layout)
          : file_name(fileName),
            format(fmt),
            layout(layout)
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
        /// Output format object which may have useful members
        FormatT                        format;
        /// Output layout
        LayoutT                        layout;
      };

      template <typename ArgT, typename FormatT, typename LayoutT>
      Sink(ArgT&&, FormatT, LayoutT) -> Sink<FormatT, LayoutT>;

      /// Variant type for all possible sink types
      using SinkVariant = std::variant<Sink<Csv, Wide>, Sink<Json, Wide>, Sink<Yaml, Wide>, Sink<Csv, Long>, Sink<Json, Long>, Sink<Yaml, Long>>;

      /**
       * @brief Factory function mapping format enum value to sink type
       */
      template <typename T>
      static SinkVariant make_sink(SinkSpec spec, T&& arg)
      {
        switch (spec.format)
        {
        case Format::CSV:
          return make_sink(std::forward<T>(arg), Csv{spec.delim}, spec.layout);
          break;
        case Format::JSON:
          return make_sink(std::forward<T>(arg), Json{}, spec.layout);
          break;
        case Format::YAML:
          return make_sink(std::forward<T>(arg), Yaml{}, spec.layout);
          break;
        }
        throw std::runtime_error("Invalid monitor output format");
      }

      /**
       * @brief Factory function mapping layout enum value to sink type
       */
      template <typename T, typename FormatT>
      static SinkVariant make_sink(T&& arg, FormatT fmt, Layout layout)
      {
        switch (layout)
        {
        case Layout::WIDE:
          return Sink(std::forward<T>(arg), fmt, Wide{});
          break;
        case Layout::LONG:
          return Sink(std::forward<T>(arg), fmt, Long{});
          break;
        }
        throw std::runtime_error("Invalid monitor output layout");
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
