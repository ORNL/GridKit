/**
 * @file ComponentDefinition.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 *
 * @brief Compile-time definition tables for component models.
 */

#pragma once

#include <array>
#include <cstddef>
#include <string_view>

#include <GridKit/Utilities/Enum.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    template <Utilities::SizedEnum EnumT>
    struct ParameterDefinition
    {
      EnumT            key;
      std::string_view units;
      std::string_view symbol;
      std::string_view description;

      constexpr std::string_view name() const
      {
        return Utilities::enum_name(key);
      }
    };

    template <Utilities::SizedEnum EnumT>
    struct PortDefinition
    {
      EnumT            key;
      std::string_view description;

      constexpr std::string_view name() const
      {
        return Utilities::enum_name(key);
      }
    };

    template <Utilities::SizedEnum EnumT>
    struct MonitorDefinition
    {
      EnumT            key;
      std::string_view units;
      std::string_view description;

      constexpr std::string_view name() const
      {
        return Utilities::enum_name(key);
      }
    };

    /// True when every enumeration value has exactly one described row
    template <typename RowT, std::size_t N>
    consteval bool rowsComplete(const std::array<RowT, N>& rows)
    {
      for (std::size_t value = 0; value < N; ++value)
      {
        std::size_t described = 0;
        for (const auto& row : rows)
        {
          if (static_cast<std::size_t>(row.key) == value && !row.description.empty())
          {
            ++described;
          }
        }
        if (described != 1)
        {
          return false;
        }
      }
      return true;
    }

    /// Per-class metadata of a component model, one table per `ComponentData` enumeration
    template <Utilities::SizedEnum ParametersT,
              Utilities::SizedEnum BusesT,
              Utilities::SizedEnum SignalInputsT,
              Utilities::SizedEnum SignalOutputsT,
              Utilities::SizedEnum MonitorableVariablesT>
    struct ComponentDefinition
    {
      template <typename RowT, typename EnumT>
      using Table = std::array<RowT, Utilities::enum_size<EnumT>()>;

      std::string_view                                                       device_class;
      std::string_view                                                       description;
      Table<ParameterDefinition<ParametersT>, ParametersT>                   parameters;
      Table<PortDefinition<BusesT>, BusesT>                                  buses;
      Table<PortDefinition<SignalInputsT>, SignalInputsT>                    signal_inputs;
      Table<PortDefinition<SignalOutputsT>, SignalOutputsT>                  signal_outputs;
      Table<MonitorDefinition<MonitorableVariablesT>, MonitorableVariablesT> monitorable_variables;

      /// True when the class and description are set and every section is complete
      consteval bool complete() const
      {
        return !device_class.empty() && !description.empty() && rowsComplete(parameters)
               && rowsComplete(buses) && rowsComplete(signal_inputs)
               && rowsComplete(signal_outputs) && rowsComplete(monitorable_variables);
      }
    };
  } // namespace PhasorDynamics
} // namespace GridKit
