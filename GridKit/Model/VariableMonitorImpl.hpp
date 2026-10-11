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

      VariableMonitor() = default;

      /**
       * @brief Construct monitor with component/bus class, identifier, and set
       * of variables selected for monitoring
       *
       * @param device_class Device class, such as `Genrou`
       * @param id Device identifier within its class
       * @param variables Variables selected for monitoring
       */
      VariableMonitor(const std::string&            device_class,
                      const std::string&            id,
                      const std::set<VariableEnum>& variables)
        : device_class_(device_class),
          id_(id)
      {
        std::ranges::copy(variables, std::back_inserter(variables_));
      }

      /**
       * @brief Construct from ObjData object
       *
       * @param data Expected to be derived from ComponentData
       */
      VariableMonitor(const ObjData& data)
        : VariableMonitor(data.device_class,
                          data.disambiguation_string,
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
        f_[static_cast<size_t>(enum_integer(v))] = ValueFormatterType<FuncT>{f};
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

        void operator()(std::string& out) const
        {
          using RetValueT = std::remove_cvref_t<RetT>;

          if constexpr (std::is_floating_point_v<RetValueT>
                        || std::is_same_v<RetValueT, ScalarT>)
          {
            VariableMonitorDetail::formatReal(out, static_cast<RealT>(f()));
          }
          else
          {
            out += (std::ostringstream{} << f()).str();
          }
        }
      };

      template <typename FuncT>
      using ValueFormatterType =
          ValueFormatterImpl<FuncT, std::invoke_result_t<FuncT>>;

      ///@}

      void appendEntries(std::vector<Entry>& entries) const override
      {
        for (auto v : variables_)
        {
          entries.emplace_back(device_class_, id_, enum_name(v), &f_[static_cast<size_t>(enum_integer(v))]);
        }
      }

    private:
      /// Compile-time constant size: length of enum value list
      static constexpr auto enum_size_ = enum_count<VariableEnum>();

      /// Set of functions associated with each enum value
      std::array<ValueFormatter, enum_size_> f_;
      /// Set of selected enum values
      std::vector<VariableEnum>              variables_;
      /// Device class, such as `Genrou`
      std::string                            device_class_;
      /// Device identifier within its class
      std::string                            id_;
    };
  } // namespace Model
} // namespace GridKit
