#pragma once

#include <concepts>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <type_traits>
#include <variant>

#include <magic_enum/magic_enum.hpp>
#include <nlohmann/json.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    using json = nlohmann::json;

    /**
     * @brief Unified interface for `Component` data containers
     *
     * @tparam RealT Real parameter data type
     * @tparam IdxT  Integer parameter data type
     */
    template <typename RealT, typename IdxT, typename Parameters, typename Ports, typename MonitorableVariables>
      requires std::is_enum_v<Parameters> && std::is_enum_v<Ports> && std::is_enum_v<MonitorableVariables>
    struct ComponentData
    {
      /// Class of device this is for
      std::string device_class;

      /// Mapping of parameters to parameter values
      std::map<Parameters, std::variant<bool, RealT, IdxT>> parameters;

      /// Mapping of ports to port values
      std::map<Ports, IdxT> ports;

      /// Set of variables being monitored
      std::set<MonitorableVariables> monitored_variables;

      std::optional<RealT> freq_base; ///< Override for the system-wide base frequency
      std::optional<RealT> va_base;   ///< Override for the system-wide power base

      std::string disambiguation_string; ///< Disambiguation string for this device

    protected:
      ComponentData() = default;
    };

    /// JSON parser function for the `ComponentData` class and descendants
    template <typename RealT,
              typename IdxT,
              typename Parameters,
              typename Ports,
              typename MonitorableVariables>
      requires std::is_enum_v<Parameters> && std::is_enum_v<Ports> && std::is_enum_v<MonitorableVariables>
    void from_json(const json& j, ComponentData<RealT, IdxT, Parameters, Ports, MonitorableVariables>& c)
    {
      j.at("class").get_to(c.device_class);

      for (auto& raw_parameter : j.at("params").items())
      {
        auto key = magic_enum::enum_cast<Parameters>(raw_parameter.key());
        if (key.has_value())
        {
          // NOTE: this is necessary because it doesn't seem like nlohmann/json handles std::variant
          //       out of the box
          if (raw_parameter.value().is_boolean())
          {
            c.parameters[key.value()] = raw_parameter.value().template get<bool>();
          }
          else if (raw_parameter.value().is_number_float())
          {
            c.parameters[key.value()] = raw_parameter.value().template get<RealT>();
          }
          else if (raw_parameter.value().is_number_integer())
          {
            c.parameters[key.value()] = raw_parameter.value().template get<IdxT>();
          }
          else
          {
            throw "Invalid initial parameter";
          }
        }
        else
        {
          throw "Invalid initial parameter";
        }
      }

      for (auto& raw_port : j.at("ports").items())
      {
        auto key = magic_enum::enum_cast<Ports>(raw_port.key());
        if (key.has_value())
        {
          raw_port.value().get_to(c.ports[key.value()]);
        }
        else
        {
          throw "Invalid port mapping";
        }
      }

      if (j.contains("freq_base"))
      {
        j.at("freq_base").get_to(c.freq_base);
      }

      if (j.contains("va_base"))
      {
        j.at("va_base").get_to(c.va_base);
      }

      j.at("id").get_to(c.disambiguation_string);

      if (j.contains("mon"))
      {
        for (auto& raw_monitored_variable : j.at("mon"))
        {
          auto monitored = magic_enum::enum_cast<MonitorableVariables>(raw_monitored_variable.get<std::string>());
          if (monitored.has_value())
          {
            c.monitored_variables.insert(monitored.value());
          }
          else
          {
            throw "Invalid monitored variable";
          }
        }
      }
    }
  } // namespace PhasorDynamics
} // namespace GridKit
