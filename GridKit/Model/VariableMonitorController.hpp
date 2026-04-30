#pragma once

/**
 * @file
 */

#include <ostream>
#include <string>

#include <GridKit/ScalarTraits.hpp>

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
    template <typename ScalarP>
    class VariableMonitorController : public VariableMonitorBase
    {
    public:
      /// Simulation scalar value type
      using ScalarT  = ScalarP;
      /// Underlying real value type
      using RealT    = typename GridKit::ScalarTraits<ScalarT>::RealT;
      ///@{
      /// @brief Alias
      using Format   = VariableMonitorFormat;
      using SinkSpec = VariableMonitorBase::SinkSpec;
      using Csv      = VariableMonitorBase::Csv;
      using Json     = VariableMonitorBase::Json;
      using Yaml     = VariableMonitorBase::Yaml;
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
          monitors_.push_back(monitor);
        }
      }

      /**
       * @brief Add output sink based on spec
       *
       * If `spec.file_name` is empty, `std::cout` is used.
       *
       * @param spec Specifies details for the sink.
       */
      void addSink(const SinkSpec& spec)
      {
        if (spec.file_name.empty())
        {
          sinks_.push_back(make_sink(spec, std::cout));
        }
        else
        {
          sinks_.push_back(make_sink(spec, spec.file_name));
        }
      }

      /**
       * @brief Provide additional top-level variables (alongside time) to be
       * printed before submonitors.
       *
       * @param label Header label for CSV; key for JSON or YAML
       * @param value Pointer to monitored variable
       */
      void addVariable(const std::string& label, const RealT* value)
      {
        variables_.emplace_back(label, value);
      }

      bool empty() const override
      {
        return (!time_) || monitors_.empty();
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

      /// @copydoc VariableMonitorBase::printHeader
      using VariableMonitorBase::printHeader;

      void printHeader(std::ostream& os, Csv csv) const override
      {
        os << "t";
        for (auto&& var : variables_)
        {
          os << csv.delim << var.label;
        }
      }

      void printHeader(std::ostream& os, Json) const override
      {
        os << "[\n";
      }

      /**
       * @brief Organize header output for this and all submonitors
       */
      template <typename FormatT>
      void printFullHeader(std::ostream& os, FormatT fmt) const
      {
        this->printHeader(os, fmt);
        for (auto* mon : monitors_)
        {
          mon->printHeader(os, fmt);
        }
        os << '\n';
      }

      /**
       * @brief Print header for all sinks
       */
      void printHeader() const
      {
        for (auto&& sink : sinks_)
        {
          std::visit([this](auto&& sink)
                     { this->printFullHeader(*sink.os, sink.format); },
                     sink);
        }
      }

      void print(std::ostream& os, Csv csv) const override
      {
        os << *time_;
        for (auto&& var : variables_)
        {
          os << csv.delim << *var.value;
        }

        for (auto* mon : monitors_)
        {
          mon->print(os, csv);
        }
      }

      void print(std::ostream& os, Json json) const override
      {
        std::string indent = "  ";
        if (json.after_first)
        {
          os << indent << ",\n";
        }
        os << indent << "{\n";
        indent.append(2, ' ');

        os << indent << std::quoted("t") << ": " << *time_ << ",\n";
        for (auto&& var : variables_)
        {
          os << indent << std::quoted(var.label) << ": " << *var.value << ",\n";
        }

        auto after_first = false;
        for (auto* mon : monitors_)
        {
          if (after_first)
          {
            os << ",\n";
          }
          mon->print(os, Json());
          after_first = true;
        }

        indent.erase(indent.size() - 2);
        os << '\n'
           << indent << "}";
      }

      void print(std::ostream& os, Yaml) const override
      {
        std::string indent = "  ";
        os << indent << "- t: " << *time_ << '\n';
        indent.append(2, ' ');
        for (auto&& var : variables_)
        {
          os << indent << var.label << ": " << *var.value << '\n';
        }

        for (auto* mon : monitors_)
        {
          mon->print(os, Yaml());
        }
      }

      /**
       * @brief Organize variable output for this and all submonitors
       */
      template <typename FormatT>
      void printFull(std::ostream& os, FormatT fmt) const
      {
        const auto     orig_prec = os.precision();
        constexpr auto max_prec  = std::numeric_limits<RealT>::digits10 + 1;
        os.precision(max_prec);
        os << std::scientific;
        this->print(os, fmt);
        os << '\n';
        os << std::defaultfloat;
        os.precision(orig_prec);
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
              this->printFull(*sink.os, sink.format);
              using T = std::remove_cvref_t<decltype(sink)>;
              if constexpr (std::is_same_v<T, Sink<Json>>)
              {
                sink.format.after_first = true;
              } },
                     sink);
        }
      }

      /// @copydoc VariableMonitorBase::printFooter
      using VariableMonitorBase::printFooter;

      void printFooter(std::ostream& os, Json) const override
      {
        os << "\n]\n";
      }

      /**
       * @brief Organize footer output for this and all submonitors
       */
      template <typename FormatT>
      void printFullFooter(std::ostream& os, FormatT format) const
      {
        for (auto* mon : monitors_)
        {
          mon->printFooter(os, format);
        }
        this->printFooter(os, format);
      }

      /**
       * @brief Print footer for all sinks
       */
      void printFooter() const
      {
        for (auto&& sink : sinks_)
        {
          std::visit([this](auto&& sink)
                     { this->printFullFooter(*sink.os, sink.format); },
                     sink);
        }
      }

    private:
      /// Time variable; printed first
      const RealT* time_{nullptr};

      /**
       * @brief Define sink for a specific output format
       */
      template <typename FormatT>
      struct Sink
      {
        Sink() = delete;

        /**
         * @brief Version for an output stream that already exists
         */
        Sink(std::ostream& out, FormatT fmt)
          : os(&out), format(fmt)
        {
        }

        /**
         * @brief Version for opening an output stream for the given file
         */
        Sink(const std::string& fileName, FormatT fmt)
          : file_name(fileName),
            format(fmt)
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
      };

      template <typename ArgT, typename FormatT>
      Sink(ArgT&&, FormatT) -> Sink<FormatT>;

      /// Variant type for all possible sink types
      using SinkVariant = std::variant<Sink<Csv>, Sink<Json>, Sink<Yaml>>;

      /**
       * @brief Factory function mapping format enum value to sink type
       */
      template <typename T>
      static SinkVariant make_sink(SinkSpec spec, T&& arg)
      {
        switch (spec.format)
        {
        case Format::CSV:
          return Sink(std::forward<T>(arg), Csv{spec.delim});
          break;
        case Format::JSON:
          return Sink(std::forward<T>(arg), Json{});
          break;
        case Format::YAML:
          return Sink(std::forward<T>(arg), Yaml{});
          break;
        }
        throw std::runtime_error("Invalid monitor output format");
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
      std::vector<SinkVariant>                sinks_;
      /// Collection of submonitors
      std::vector<const VariableMonitorBase*> monitors_;

      /**
       * @brief Key/Value object for extra top-level variables
       */
      struct Variable
      {
        /// Header label for CSV; key for JSON or YAML
        std::string  label;
        /// Pointer to monitored variable
        const RealT* value;
      };

      /// Collection of extra top-level monitored variables
      std::vector<Variable> variables_;
    };
  } // namespace Model
} // namespace GridKit
