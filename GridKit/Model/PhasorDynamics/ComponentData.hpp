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
              typename Ports,
              typename MonitorableVariables>
      requires std::is_enum_v<Parameters>
               && std::is_enum_v<Ports>
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

      /// Mapping of ports to port values
      std::map<Ports, IdxT> ports;

      /// Mapping of signal input ports to fixed input delays.
      std::map<Ports, RealT> port_delays;

      /// Set of variables being monitored
      std::set<MonitorableVariables> monitored_variables;

      std::string disambiguation_string; ///< Disambiguation string for this device

    protected:
      ComponentData() = default;
    };
  } // namespace PhasorDynamics
} // namespace GridKit
