#pragma once

#include <magic_enum/magic_enum.hpp>

#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit
{
  namespace Model
  {
    template <typename ScalarT>
    class VariableMonitorController;

    using magic_enum::enum_count;
    using magic_enum::enum_integer;
    using magic_enum::enum_name;

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
      template <typename>
      friend class VariableMonitorController;

    public:
      /// Underlying real value type
      using RealT        = typename GridKit::ScalarTraits<ScalarT>::RealT;
      /// Type of (EvalT)Data class expected to have MonitorableVariables enum
      using ObjData      = DataT<RealT, IdxT>;
      /// Enum of valid monitorable variables
      using VariableEnum = typename ObjData::MonitorableVariables;
      ///@{
      /// @brief Alias
      using Csv          = VariableMonitorBase::Csv;
      using Json         = VariableMonitorBase::Json;
      using Yaml         = VariableMonitorBase::Yaml;
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
      template <typename FuncT>
      void set(VariableEnum v, FuncT f)
      {
        f_[static_cast<size_t>(enum_integer(v))] = ValuePrinter{f};
      }

    private:
      ///@{
      /**
       * @brief Functors to handle printing different types
       */
      template <typename FuncT, typename RetT>
      struct ValuePrinterImpl
      {
        FuncT f;

        void operator()(std::ostream& os) const
        {
          os << f();
        }

        void appendCsv(std::string& out) const
        {
          using RetValueT = std::remove_cvref_t<RetT>;

          if constexpr (std::is_floating_point_v<RetValueT>)
          {
            VariableMonitorDetail::appendCsvReal(out, static_cast<RealT>(f()));
          }
          else
          {
            std::ostringstream os;
            os << f();
            out += os.str();
          }
        }
      };

      template <typename FuncT>
      struct ValuePrinterImpl<FuncT, ScalarT>
      {
        FuncT f;

        void operator()(std::ostream& os) const
        {
          os << static_cast<RealT>(f());
        }

        void appendCsv(std::string& out) const
        {
          VariableMonitorDetail::appendCsvReal(out, static_cast<RealT>(f()));
        }
      };

      template <typename FuncT>
      using ValuePrinterType =
          ValuePrinterImpl<FuncT, std::invoke_result_t<FuncT>>;

      class ValuePrinter
      {
      public:
        ValuePrinter() = default;

        template <typename FuncT>
        ValuePrinter(FuncT f)
        {
          auto printer = ValuePrinterType<FuncT>{f};
          impl_        = [printer](std::ostream& os)
          {
            printer(os);
          };
          append_csv_ = [printer](std::string& out)
          {
            printer.appendCsv(out);
          };
        }

        void appendCsv(std::string& out) const
        {
          append_csv_(out);
        }

      private:
        std::function<void(std::ostream&)> impl_;
        std::function<void(std::string&)>  append_csv_;

        friend std::ostream& operator<<(std::ostream& os, const ValuePrinter& p)
        {
          p.impl_(os);
          return os;
        }
      };

      ///@}

      void printHeader(std::ostream& os, Csv csv) const override
      {
        for (auto v : variables_)
        {
          os << csv.delim << label_ << '_' << enum_name(v);
        }
      }

      void print(std::ostream& os, Csv csv) const override
      {
        for (auto v : variables_)
        {
          os << csv.delim << f_[static_cast<size_t>(enum_integer(v))];
        }
      }

      void appendCsvHeader(std::string& out, Csv csv) const override
      {
        for (auto v : variables_)
        {
          out += csv.delim;
          out += label_;
          out += '_';
          out += enum_name(v);
        }
      }

      void appendCsv(std::string& out, Csv csv) const override
      {
        for (auto v : variables_)
        {
          out += csv.delim;
          f_[static_cast<size_t>(enum_integer(v))].appendCsv(out);
        }
      }

      /**
       * @brief Print single variable
       */
      void print(std::ostream& os, VariableEnum v, Json) const
      {
        os << indent_ << std::quoted(enum_name(v)) << ": "
           << f_[static_cast<size_t>(enum_integer(v))] << ",\n";
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
        os << indent_ << enum_name(v) << ": " << f_[static_cast<size_t>(enum_integer(v))]
           << '\n';
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
      static constexpr auto enum_size_ = enum_count<VariableEnum>();

      /// Set of functions associated with each enum value
      std::array<ValuePrinter, enum_size_> f_;
      /// Set of selected enum values
      std::vector<VariableEnum>            variables_;
      /// Indent string used for formatting
      mutable std::string                  indent_{"    "};
      /// Monitor disambiguation label
      std::string                          label_;
    };
  } // namespace Model
} // namespace GridKit
