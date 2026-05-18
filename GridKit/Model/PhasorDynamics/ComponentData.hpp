#pragma once

#include <cstddef>
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
    /// Default external-variable enum for components without device init support.
    enum class NoExternalVariables : std::size_t
    {
      MAXIMUM
    };

    /**
     * @brief Unified interface for `Component` data containers
     *
     * @tparam RealT Real parameter data type
     * @tparam IdxT  Integer parameter data type
     */
    template <typename RealT,
              typename IdxT,
              typename Parameters,
              typename Ports,
              typename MonitorableVariables,
              typename ExternalVariables = NoExternalVariables>
      requires std::is_enum_v<Parameters>
               && std::is_enum_v<Ports>
               && std::is_enum_v<MonitorableVariables>
               && std::is_enum_v<ExternalVariables>
               && requires { ExternalVariables::MAXIMUM; }
    struct ComponentData
    {
      using ValueT = std::variant<bool, RealT, IdxT>;

      /// Class of device this is for
      std::string device_class;

      /// Mapping of parameters to parameter values
      std::map<Parameters, ValueT> parameters;

      /// Mapping of external variables to initialization values
      std::map<ExternalVariables, ValueT> external_initial_values;

      /// Mapping of ports to port values
      std::map<Ports, IdxT> ports;

      /// Set of variables being monitored
      std::set<MonitorableVariables> monitored_variables;

      /// Number of invalid input records detected while parsing this component
      std::size_t input_error_count{0};

      std::string disambiguation_string; ///< Disambiguation string for this device

    protected:
      ComponentData() = default;
    };
  } // namespace PhasorDynamics
} // namespace GridKit
