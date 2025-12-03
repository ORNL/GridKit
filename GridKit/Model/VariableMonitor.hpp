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

    template <typename ScalarT, typename IdxT, template <typename, typename> typename EvalT, template <typename, typename> typename DataT>
    class VariableMonitor<EvalT<ScalarT, IdxT>, DataT> : public VariableMonitorBase
    {
    public:
      using RealT         = typename GridKit::ScalarTraits<ScalarT>::RealT;
      using ObjData       = DataT<RealT, IdxT>;
      using VariableEnum  = typename ObjData::MonitorableVariables;
      using EnumId        = std::underlying_type_t<VariableEnum>;
      using ValueFunction = std::function<ScalarT(void)>;
      using Csv           = VariableMonitorBase::Csv;
      using Json          = VariableMonitorBase::Json;
      using Yaml          = VariableMonitorBase::Yaml;

      VariableMonitor() = default;

      VariableMonitor(const std::string&            label,
                      const std::set<VariableEnum>& variables)
        : label_(label)
      {
        std::ranges::copy(variables, std::back_inserter(variables_));
      }

      VariableMonitor(const ObjData& data)
        : VariableMonitor(data.device_class + "_" + data.disambiguation_string,
                          data.monitored_variables)
      {
      }

      virtual ~VariableMonitor()
      {
      }

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

      void print(std::ostream& os, VariableEnum v, Yaml) const
      {
        os << indent_ << indent_ << enumLabel(v) << ": " << f(v) << '\n';
      }

      void print(std::ostream& os, Yaml) const override
      {
        if (empty())
        {
          return;
        }
        os << indent_ << label_ << ":\n";
        for (auto v : variables_)
        {
          print(os, v, Yaml());
        }
      }

      bool empty() const override
      {
        return variables_.empty();
      }

      void set(VariableEnum v, ValueFunction f)
      {
        f_[Utilities::enumId(v)] = f;
      }

    private:
      auto f(VariableEnum v) const
      {
        return static_cast<RealT>(f_[Utilities::enumId(v)]());
      }

      static constexpr auto                 enum_size_ = Utilities::enumSize<VariableEnum>;
      std::array<ValueFunction, enum_size_> f_;
      std::vector<VariableEnum>             variables_;
      mutable std::string                   indent_{"    "};
      std::string                           label_;
    };

    template <typename ScalarT, typename IdxT>
    class VariableMonitor<PhasorDynamics::SystemModel<ScalarT, IdxT>, PhasorDynamics::SystemModelData> : public VariableMonitorBase
    {
    public:
      using RealT         = typename GridKit::ScalarTraits<ScalarT>::RealT;
      using ValueFunction = std::function<RealT(void)>;
      using Format        = VariableMonitorFormat;
      using SinkSpec      = VariableMonitorBase::SinkSpec;
      using Csv           = VariableMonitorBase::Csv;
      using Json          = VariableMonitorBase::Json;
      using Yaml          = VariableMonitorBase::Yaml;

      VariableMonitor() = default;

      explicit VariableMonitor(const RealT& time_var)
        : time_(&time_var)
      {
      }

      virtual ~VariableMonitor()
      {
      }

      void addMonitor(const VariableMonitorBase* monitor)
      {
        if (monitor && !monitor->empty())
        {
          monitors_.push_back(monitor);
        }
      }

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

      void addVariable(const std::string& label, const RealT* value)
      {
        variables_.emplace_back(label, value);
      }

      bool empty() const
      {
        return (!time_) || monitors_.empty();
      }

      void start() const
      {
        if (!empty())
        {
          printHeader();
        }
      }

      void stop() const
      {
        if (!empty())
        {
          printFooter();
        }
      }

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

      using VariableMonitorBase::printFooter;

      void printFooter(std::ostream& os, Json) const override
      {
        os << "\n]\n";
      }

      template <typename FormatT>
      void printFullFooter(std::ostream& os, FormatT format) const
      {
        for (auto* mon : monitors_)
        {
          mon->printFooter(os, format);
        }
        this->printFooter(os, format);
      }

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
      const RealT* time_{nullptr};

      template <typename FormatT>
      struct Sink
      {
        Sink() = delete;

        Sink(std::ostream& out, FormatT fmt)
          : os(out), format(fmt)
        {
        }

        Sink(const std::string& fileName, FormatT fmt)
          : file_stream(std::make_unique<std::ofstream>(fileName)),
            os(*file_stream),
            format(fmt)
        {
        }

        Sink(Sink&&) = default;

        std::unique_ptr<std::ofstream> file_stream;
        std::ostream&                  os;
        FormatT                        format;
      };

      using SinkVariant = std::variant<Sink<Csv>, Sink<Json>, Sink<Yaml>>;

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

      std::vector<SinkVariant>                sinks_;
      std::vector<const VariableMonitorBase*> monitors_;

      struct Variable
      {
        std::string  label;
        const RealT* value;
      };

      std::vector<Variable> variables_;

      mutable bool after_first_{false};
    };

  } // namespace Model
} // namespace GridKit
