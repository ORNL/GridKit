#pragma once

#include <algorithm>
#include <array>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <set>
#include <string>
#include <type_traits>
#include <variant>
#include <vector>

#include <GridKit/Model/VariableMonitorBase.hpp>
#include <GridKit/ScalarTraits.hpp>
#include <GridKit/Utilities/Enum.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    template <typename, typename>
    class SystemModel;

    template <typename, typename>
    class SystemModelData;
  } // namespace PhasorDynamics

  namespace Model
  {
    template <typename EvalT, template <typename, typename> typename DataT>
    class VariableMonitor
    {
    };

    /**
     * @brief Manage printing variables associated with a specific component or
     * bus
     *
     * Implementations of the print functions print under the assumption of a
     * larger context. For example, the Csv version prints the delimiter and the
     * value (or label for the header) for each monitored value without a line
     * break. That way other monitors can print likewise on the same line, and
     * the line can be ended by the control monitor.
     */
    template <typename ScalarT,
              typename IdxT,
              template <typename, typename> typename EvalT,
              template <typename, typename> typename DataT>
    class VariableMonitor<EvalT<ScalarT, IdxT>, DataT>
      : public VariableMonitorBase
    {
      template <typename, template <typename, typename> typename>
      friend class VariableMonitor;

    public:
      /// Underlying real value type
      using RealT         = typename GridKit::ScalarTraits<ScalarT>::RealT;
      /// Type of (EvalT)Data class expected to have MonitorableVariables enum
      using ObjData       = DataT<RealT, IdxT>;
      /// Enum of valid monitorable variables
      using VariableEnum  = typename ObjData::MonitorableVariables;
      /// Abstraction type for functions returning a monitored value
      using ValueFunction = std::function<ScalarT(void)>;
      ///@{
      /// @brief Alias
      using Csv           = VariableMonitorBase::Csv;
      using Json          = VariableMonitorBase::Json;
      using Yaml          = VariableMonitorBase::Yaml;
      ///@}

      VariableMonitor() = default;

      /**
       * @brief Construct monitor with component/bus label and set of variables
       * selected for monitoring
       *
       * @param label Monitor label used to disambiguate like variables from
       * different objects of same type
       * @param variables Variables selected for monitoring
       */
      VariableMonitor(const std::string&            label,
                      const std::set<VariableEnum>& variables)
        : label_(label)
      {
        std::ranges::copy(variables, std::back_inserter(variables_));
      }

      /**
       * @brief Construct from ObjData object
       *
       * Constructs the monitor label from elements of the data object
       *
       * @param data Expected to be derived from ComponentData
       */
      VariableMonitor(const ObjData& data)
        : VariableMonitor(data.device_class + "_" + data.disambiguation_string,
                          data.monitored_variables)
      {
      }

      virtual ~VariableMonitor()
      {
      }

      bool empty() const override
      {
        return variables_.empty();
      }

      /**
       * @brief Associate a value getter with a variable enum value
       *
       * @note This does not designate the variable for printing. It defines how
       * to get the variable if it is printed.
       */
      void set(VariableEnum v, ValueFunction f)
      {
        f_[Utilities::enumId(v)] = f;
      }

    protected:
      void printHeader(std::ostream& os, Csv csv) const override
      {
        for (auto v : variables_)
        {
          os << csv.delim << label_ << '_' << enumLabel(v);
        }
      }

      void print(std::ostream& os, Csv csv) const override
      {
        for (auto v : variables_)
        {
          os << csv.delim << f(v);
        }
      }

      /**
       * @brief Print single variable
       */
      void print(std::ostream& os, VariableEnum v, Json) const
      {
        os << indent_ << std::quoted(enumLabel(v)) << ": " << f(v) << ",\n";
      }

      void print(std::ostream& os, Json) const override
      {
        if (empty())
        {
          return;
        }
        os << indent_ << std::quoted(label_) << ": {\n";
        indent_.append(2, ' ');
        std::ostringstream v_os;
        v_os.copyfmt(os);
        for (auto v : variables_)
        {
          print(v_os, v, Json());
        }
        auto vars = v_os.view();
        vars.remove_suffix(2);
        os << vars << '\n';
        indent_.erase(indent_.size() - 2);
        os << indent_ << "}";
      }

      void printHeader(std::ostream&, Yaml) const override
      {
      }

      /**
       * @brief Print single variable
       */
      void print(std::ostream& os, VariableEnum v, Yaml) const
      {
        os << indent_ << enumLabel(v) << ": " << f(v) << '\n';
      }

      void print(std::ostream& os, Yaml) const override
      {
        if (empty())
        {
          return;
        }
        os << indent_ << label_ << ":\n";
        indent_.append(2, ' ');
        for (auto v : variables_)
        {
          print(os, v, Yaml());
        }
        indent_.erase(indent_.size() - 2);
      }

    private:
      /// Compile-time constant size: length of enum value list
      static constexpr auto enum_size_ = Utilities::enumSize<VariableEnum>;

      /**
       * @brief Convenience function to access value associated with enum value
       */
      auto f(VariableEnum v) const
      {
        return static_cast<RealT>(f_[Utilities::enumId(v)]());
      }

      /// Set of functions associated with each enum value
      std::array<ValueFunction, enum_size_> f_;
      /// Set of selected enum values
      std::vector<VariableEnum>             variables_;
      /// Indent string used for formatting
      mutable std::string                   indent_{"    "};
      /// Monitor disambiguation label
      std::string                           label_;
    };

    /**
     * @brief Monitor associated with the SystemModel; controls output from
     * component and bus monitors.
     *
     * High-level print functions (without parameters) manage printing for all
     * monitors for multiple output sinks.
     */
    template <typename ScalarT, typename IdxT>
    class VariableMonitor<PhasorDynamics::SystemModel<ScalarT, IdxT>,
                          PhasorDynamics::SystemModelData>
      : public VariableMonitorBase
    {
    public:
      /// Underlying real value type
      using RealT         = typename GridKit::ScalarTraits<ScalarT>::RealT;
      /// Abstraction type for functions returning a monitored value
      using ValueFunction = std::function<RealT(void)>;
      ///@{
      /// @brief Alias
      using Format        = VariableMonitorFormat;
      using SinkSpec      = VariableMonitorBase::SinkSpec;
      using Csv           = VariableMonitorBase::Csv;
      using Json          = VariableMonitorBase::Json;
      using Yaml          = VariableMonitorBase::Yaml;
      ///@}

      /// Default to empty monitor
      VariableMonitor() = default;

      /**
       * @brief Constructor expects a time variable to monitor
       */
      explicit VariableMonitor(const RealT& time_var)
        : time_(&time_var)
      {
      }

      virtual ~VariableMonitor()
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
      void start() const
      {
        if (!empty())
        {
          printHeader();
        }
      }

      /**
       * @brief Print footer if we're monitoring
       */
      void stop() const
      {
        if (!empty())
        {
          printFooter();
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
                     { printFullHeader(sink.os, sink.format); },
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
              printFull(sink.os, sink.format);
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
                     { printFullFooter(sink.os, sink.format); },
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
          : os(out), format(fmt)
        {
        }

        /**
         * @brief Version for opening an output stream for the given file
         */
        Sink(const std::string& fileName, FormatT fmt)
          : file_stream(std::make_unique<std::ofstream>(fileName)),
            os(*file_stream),
            format(fmt)
        {
        }

        /**
         * @brief Have to move because of unique_ptr
         */
        Sink(Sink&&) = default;

        /// Output file stream (if we opened one)
        std::unique_ptr<std::ofstream> file_stream;
        /// Output stream for printing
        std::ostream&                  os;
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
