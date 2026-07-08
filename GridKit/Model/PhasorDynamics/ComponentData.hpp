#pragma once

#include <map>
#include <optional>
#include <set>
#include <string>
#include <type_traits>
#include <variant>

namespace GridKit
{
  namespace PhasorDynamics
  {
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
      using IdxT  = index_type;

      using Parameters           = parameters_enum;
      using Buses                = buses_enum;
      using SignalInputs         = signal_inputs_enum;
      using SignalOutputs        = signal_outputs_enum;
      using MonitorableVariables = monitorable_variables_enum;

      /// Class of device this is for
      std::string device_class;

      /// Mapping of parameters to parameter values
      std::map<Parameters, std::variant<bool, RealT, IdxT>> parameters;

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
