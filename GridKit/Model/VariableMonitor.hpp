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
#include <vector>

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
    enum class VariableMonitorFormat
    {
      CSV,
      JSON,
      YAML
    };

    class VariableMonitorBase
    {
    public:
      struct Csv
      {
      };

      struct Json
      {
      };

      struct Yaml
      {
      };

      using Format = VariableMonitorFormat;

      struct SinkSpec
      {
        std::string file_name;
        Format      format;
      };

      virtual ~VariableMonitorBase()
      {
      }

      virtual bool empty() const = 0;

      virtual void printHeader(std::ostream&, Csv) const = 0;
      virtual void print(std::ostream&, Csv) const       = 0;

      virtual void printFooter(std::ostream&, Csv) const
      {
      }

      virtual void printHeader(std::ostream&, Json) const
      {
      }

      virtual void print(std::ostream&, Json) const = 0;

      virtual void printFooter(std::ostream&, Json) const
      {
      }

      virtual void printHeader(std::ostream&, Yaml) const
      {
      }

      virtual void print(std::ostream&, Yaml) const = 0;

      virtual void printFooter(std::ostream&, Yaml) const
      {
      }
    };

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

      void printHeader(std::ostream& os, Csv) const override
      {
        for (auto v : variables_)
        {
          os << delim_ << label_ << '_' << enumLabel(v);
        }
      }

      void print(std::ostream& os, Csv) const override
      {
        for (auto v : variables_)
        {
          os << delim_ << f(v);
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
      std::string                           delim_{","};
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
        sinks_.emplace_back(std::cout, Format::CSV);
      }

      virtual ~VariableMonitor()
      {
        stop();
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
          sinks_.front().format = spec.format;
        }
        else
        {
          sinks_.emplace_back(spec.file_name, spec.format);
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

      void printHeader(std::ostream& os, Csv) const override
      {
        os << "t";
        for (auto&& var : variables_)
        {
          os << delim_ << var.label;
        }
      }

      void printHeader(std::ostream& os, Json) const override
      {
        os << "[\n";
      }

      template <typename FormatT>
      void printFullHeader(std::ostream& os, FormatT format) const
      {
        this->printHeader(os, format);
        for (auto* mon : monitors_)
        {
          mon->printHeader(os, format);
        }
        os << '\n';
      }

      void printHeader() const
      {
        for (auto&& sink : sinks_)
        {
          switch (sink.format)
          {
          case Format::CSV:
            printFullHeader(sink.os, Csv());
            break;
          case Format::JSON:
            printFullHeader(sink.os, Json());
            break;
          case Format::YAML:
            printFullHeader(sink.os, Yaml());
            break;
          }
        }
      }

      void print(std::ostream& os, Csv) const override
      {
        os << *time_;
        for (auto&& var : variables_)
        {
          os << delim_ << *var.value;
        }

        for (auto* mon : monitors_)
        {
          mon->print(os, Csv());
        }
      }

      void print(std::ostream& os, Json) const override
      {
        static bool after_first = false;

        std::string indent = "  ";
        if (after_first)
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

        after_first = false;
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

        after_first = true;
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
      void printFull(std::ostream& os, FormatT format) const
      {
        const auto     orig_prec = os.precision();
        constexpr auto max_prec  = std::numeric_limits<RealT>::digits10 + 1;
        os.precision(max_prec);
        os << std::scientific;
        this->print(os, format);
        os << '\n';
        os << std::defaultfloat;
        os.precision(orig_prec);
      }

      void print() const
      {
        for (auto&& sink : sinks_)
        {
          switch (sink.format)
          {
          case Format::CSV:
            printFull(sink.os, Csv());
            break;
          case Format::JSON:
            printFull(sink.os, Json());
            break;
          case Format::YAML:
            printFull(sink.os, Yaml());
            break;
          }
        }
      }

      using VariableMonitorBase::printFooter;

      void printFooter(std::ostream& os, Json) const override
      {
        os << "]";
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
          switch (sink.format)
          {
          case Format::CSV:
            printFullFooter(sink.os, Csv());
            break;
          case Format::JSON:
            printFullFooter(sink.os, Json());
            break;
          case Format::YAML:
            printFullFooter(sink.os, Yaml());
            break;
          }
        }
      }

    private:
      const RealT* time_{nullptr};

      struct Sink
      {
        Sink() = delete;

        Sink(std::ostream& out, Format fmt)
          : os(out), format(fmt)
        {
        }

        Sink(const std::string& fileName, Format fmt)
          : file_stream(std::make_unique<std::ofstream>(fileName)),
            os(*file_stream),
            format(fmt)
        {
        }

        std::unique_ptr<std::ofstream> file_stream;
        std::ostream&                  os;
        Format                         format;
      };

      std::vector<Sink>                       sinks_;
      std::vector<const VariableMonitorBase*> monitors_;

      struct Variable
      {
        std::string  label;
        const RealT* value;
      };

      std::vector<Variable> variables_;

      std::string delim_{","};
    };

  } // namespace Model
} // namespace GridKit
