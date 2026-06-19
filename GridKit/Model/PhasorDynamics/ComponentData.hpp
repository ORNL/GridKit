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
              typename parameters_type,
              typename buses_type,
              typename signal_inputs_type,
              typename signal_outputs_type,
              typename monitorable_variables_type>
      requires std::is_enum_v<parameters_type>
               && std::is_enum_v<buses_type>
               && std::is_enum_v<signal_inputs_type>
               && std::is_enum_v<signal_outputs_type>
               && std::is_enum_v<monitorable_variables_type>
    struct ComponentData
    {
      /// Real value type
      using RealT = real_type;
      /// Index type
      using IdxT  = index_type;

      /// Parameters enum
      using Parameters           = parameters_type;
      /// Buses enum
      using Buses                = buses_type;
      /// Signal inputs enum
      using SignalInputs         = signal_inputs_type;
      /// Signal outputs enum
      using SignalOutputs        = signal_outputs_type;
      /// Monitorable variables enum
      using MonitorableVariables = monitorable_variables_type;

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

    protected:
      ComponentData() = default;
    };
  } // namespace PhasorDynamics
} // namespace GridKit
