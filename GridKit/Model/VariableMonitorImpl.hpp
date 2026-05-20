#pragma once

#include <magic_enum/magic_enum.hpp>

#include <GridKit/Model/VariableMonitor.hpp>

namespace GridKit
{
  namespace Model
  {
    template <typename scalar_type>
    class VariableMonitorController;

    using magic_enum::enum_count;
    using magic_enum::enum_integer;
    using magic_enum::enum_name;

    /**
     * @brief Manage printing variables associated with a specific component or
     * bus
     *
     * Append functions assume a larger output context. For example, the CSV
     * version appends the delimiter and the value (or label for the header) for
     * each monitored value without a line break. That way other monitors can
     * append likewise on the same line, and the line can be ended by the control
     * monitor.
     */
    template <typename scalar_type,
              typename index_type,
              template <typename, typename> typename eval_type,
              template <typename, typename> typename model_data_type>
    class VariableMonitor<eval_type<scalar_type, index_type>, model_data_type>
      : public VariableMonitorBase
    {
      template <typename>
      friend class VariableMonitorController;

    public:
      /// Underlying scalar value type
      using ScalarT      = scalar_type;
      /// Index type
      using IdxT         = index_type;
      /// Underlying real value type
      using RealT        = typename GridKit::ScalarTraits<ScalarT>::RealT;
      /// Type of (EvalT)Data class expected to have MonitorableVariables enum
      using ObjData      = model_data_type<RealT, IdxT>;
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
        f_[static_cast<size_t>(enum_integer(v))] = ValueFormatter{f};
      }

    private:
      ///@{
      /**
       * @brief Functors to handle printing different types
       */
      template <typename FuncT, typename RetT>
      struct ValueFormatterImpl
      {
        FuncT f;

        std::string operator()() const
        {
          using RetValueT = std::remove_cvref_t<RetT>;

          if constexpr (std::is_floating_point_v<RetValueT>
                        || std::is_same_v<RetValueT, ScalarT>)
          {
            return VariableMonitorDetail::formatReal(static_cast<RealT>(f()));
          }
          else
          {
            return (std::ostringstream{} << f()).str();
          }
        }
      };

      template <typename FuncT>
      using ValueFormatterType =
          ValueFormatterImpl<FuncT, std::invoke_result_t<FuncT>>;

      class ValueFormatter
      {
      public:
        ValueFormatter() = default;

        template <typename FuncT>
        ValueFormatter(FuncT f)
          : impl_{ValueFormatterType<FuncT>{f}}
        {
        }

        std::string operator()() const
        {
          return impl_();
        }

      private:
        std::function<std::string(void)> impl_;
      };

      ///@}

      void appendHeader(std::string& out, Csv csv) const override
      {
        for (auto v : variables_)
        {
          out += csv.delim;
          out += label_;
          out += '_';
          out.append(enum_name(v));
        }
      }

      void append(std::string& out, Csv csv) const override
      {
        for (auto v : variables_)
        {
          out += csv.delim;
          out += f_[static_cast<size_t>(enum_integer(v))]();
        }
      }

      /**
       * @brief Append single variable
       */
      void appendVariable(std::string& out, VariableEnum v, Json) const
      {
        out += indent_ + "\"";
        out.append(enum_name(v));
        out += "\": ";
        out += f_[static_cast<size_t>(enum_integer(v))]();
        out += ",\n";
      }

      void append(std::string& out, Json) const override
      {
        if (empty())
        {
          return;
        }

        out += indent_ + "\"" + label_ + "\": {\n";
        indent_.append(2, ' ');
        for (auto v : variables_)
        {
          appendVariable(out, v, Json());
        }
        // Remove comma after last entry
        out.erase(out.size() - 2);
        out += '\n';
        indent_.erase(indent_.size() - 2);
        out += indent_ + "}";
      }

      /**
       * @brief Append single variable
       */
      void appendVariable(std::string& out, VariableEnum v, Yaml) const
      {
        out += indent_;
        out.append(enum_name(v));
        out += ": ";
        out += f_[static_cast<size_t>(enum_integer(v))]();
        out += '\n';
      }

      void append(std::string& out, Yaml) const override
      {
        if (empty())
        {
          return;
        }

        out += indent_ + label_ + ":\n";
        indent_.append(2, ' ');
        for (auto v : variables_)
        {
          appendVariable(out, v, Yaml());
        }
        indent_.erase(indent_.size() - 2);
      }

    private:
      /// Compile-time constant size: length of enum value list
      static constexpr auto enum_size_ = enum_count<VariableEnum>();

      /// Set of functions associated with each enum value
      std::array<ValueFormatter, enum_size_> f_;
      /// Set of selected enum values
      std::vector<VariableEnum>              variables_;
      /// Indent string used for formatting
      mutable std::string                    indent_{"    "};
      /// Monitor disambiguation label
      std::string                            label_;
    };
  } // namespace Model
} // namespace GridKit
