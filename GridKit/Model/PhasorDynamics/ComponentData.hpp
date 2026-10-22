#pragma once

#include <array>
#include <cstddef>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <string_view>
#include <type_traits>
#include <variant>

#include <GridKit/Utilities/Enum.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    /// Manifest entry describing one parameter
    template <Utilities::SizedEnum EnumT>
    struct ParameterEntry
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

    /// Manifest entry describing one bus or signal port
    template <Utilities::SizedEnum EnumT>
    struct PortEntry
    {
      EnumT            key;
      std::string_view description;

      constexpr std::string_view name() const
      {
        return Utilities::enum_name(key);
      }
    };

    /// Manifest entry describing one monitored variable
    template <Utilities::SizedEnum EnumT>
    struct MonitorEntry
    {
      EnumT            key;
      std::string_view units;
      std::string_view description;

      constexpr std::string_view name() const
      {
        return Utilities::enum_name(key);
      }
    };

    /// True when every enumeration value has exactly one described entry
    template <typename EntryT, std::size_t N>
    consteval bool entriesComplete(const std::array<EntryT, N>& entries)
    {
      for (std::size_t value = 0; value < N; ++value)
      {
        std::size_t described = 0;
        for (const auto& entry : entries)
        {
          if (static_cast<std::size_t>(entry.key) == value && !entry.description.empty())
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

    /// Units, symbols and descriptions of a component class, one entry per enumeration value
    template <Utilities::SizedEnum ParametersT,
              Utilities::SizedEnum BusesT,
              Utilities::SizedEnum InputsT,
              Utilities::SizedEnum OutputsT,
              Utilities::SizedEnum MonitorsT>
    struct ComponentManifest
    {
      template <typename EntryT, typename EnumT>
      using Table = std::array<EntryT, Utilities::enum_size<EnumT>()>;

      std::string_view                                class_name;
      std::string_view                                description;
      Table<ParameterEntry<ParametersT>, ParametersT> parameters;
      Table<PortEntry<BusesT>, BusesT>                buses;
      Table<PortEntry<InputsT>, InputsT>              inputs;
      Table<PortEntry<OutputsT>, OutputsT>            outputs;
      Table<MonitorEntry<MonitorsT>, MonitorsT>       monitors;

      /// True when the class and description are set and every table is complete
      consteval bool complete() const
      {
        return !class_name.empty() && !description.empty() && entriesComplete(parameters)
               && entriesComplete(buses) && entriesComplete(inputs)
               && entriesComplete(outputs) && entriesComplete(monitors);
      }
    };

    /**
     * @brief Unified interface for `Component` data containers
     *
     * @tparam real_type  Real parameter data type
     * @tparam index_type Integer parameter data type
     */
    template <typename real_type,
              typename index_type,
              typename parameters_enum,
              typename buses_enum,
              typename signal_inputs_enum,
              typename signal_outputs_enum,
              typename monitorable_variables_enum>
      requires std::is_enum_v<parameters_enum>
               && std::is_enum_v<buses_enum>
               && std::is_enum_v<signal_inputs_enum>
               && std::is_enum_v<signal_outputs_enum>
               && std::is_enum_v<monitorable_variables_enum>
    struct ComponentData
    {
      /// Real value type
      using RealT = real_type;

      /// Index type
      using IdxT = index_type;

      /// Enumeration over parameters.
      using Parameters = parameters_enum;

      /// Enumeration over buses.
      using Buses = buses_enum;

      /// Enumeration over signal inputs.
      using SignalInputs = signal_inputs_enum;

      /// Enumeration over signal outputs.
      using SignalOutputs = signal_outputs_enum;

      /// Enumeration over monitorable variables.
      using MonitorableVariables = monitorable_variables_enum;

      /// Units, symbols and descriptions of this component class
      using Manifest = ComponentManifest<Parameters,
                                         Buses,
                                         SignalInputs,
                                         SignalOutputs,
                                         MonitorableVariables>;

      /// Class of device this is for
      std::string device_class;

      /// Mapping of parameters to parameter values
      std::map<Parameters, std::variant<std::string, bool, RealT, IdxT>> parameters;

      /// Mapping of terminal attachments to bus identifiers
      std::map<Buses, IdxT> buses;

      /// Mapping of signal inputs to signal identifiers
      std::map<SignalInputs, IdxT> signal_inputs;

      /// Mapping of signal outputs to signal identifiers
      std::map<SignalOutputs, IdxT> signal_outputs;

      /// Set of variables being monitored
      std::set<MonitorableVariables> monitored_variables;

      /// Disambiguation string for this device
      std::string disambiguation_string;
    };
  } // namespace PhasorDynamics
} // namespace GridKit
