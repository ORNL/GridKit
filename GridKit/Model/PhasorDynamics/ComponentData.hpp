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
              typename Parameters,
              typename Terminals,
              typename InputPorts,
              typename OutputPorts,
              typename MonitorableVariables>
      requires std::is_enum_v<Parameters>
               && std::is_enum_v<Terminals>
               && std::is_enum_v<InputPorts>
               && std::is_enum_v<OutputPorts>
               && std::is_enum_v<MonitorableVariables>
    struct ComponentData
    {
      /// Real value type
      using RealT = real_type;
      /// Index type
      using IdxT  = index_type;

      /// Class of device this is for
      std::string device_class;

      /// Mapping of parameters to parameter values
      std::map<Parameters, std::variant<bool, RealT, IdxT>> parameters;

      /// Mapping of terminal attachments to bus identifiers
      std::map<Terminals, IdxT> terminals;

      /// Mapping of signal input ports to signal identifiers
      std::map<InputPorts, IdxT> input_ports;

      /// Mapping of signal output ports to signal identifiers
      std::map<OutputPorts, IdxT> output_ports;

      /// Set of variables being monitored
      std::set<MonitorableVariables> monitored_variables;

      std::string disambiguation_string; ///< Disambiguation string for this device

    protected:
      ComponentData() = default;
    };
  } // namespace PhasorDynamics
} // namespace GridKit
