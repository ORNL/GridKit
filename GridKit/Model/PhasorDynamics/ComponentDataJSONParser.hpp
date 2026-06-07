#pragma once

#include <sstream>
#include <stdexcept>
#include <string>

#include <magic_enum/magic_enum.hpp>
#include <nlohmann/json.hpp>

#include <GridKit/Model/PhasorDynamics/ComponentData.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    using json = nlohmann::json;
    using Log  = ::GridKit::Utilities::Logger;

    /// JSON parser function for the `ComponentData` class and descendants
    template <typename RealT,
              typename IdxT,
              typename Parameters,
              typename Ports,
              typename MonitorableVariables>
      requires std::is_enum_v<Parameters>
               && std::is_enum_v<Ports>
               && std::is_enum_v<MonitorableVariables>
    void from_json(const json& j, ComponentData<RealT, IdxT, Parameters, Ports, MonitorableVariables>& c)
    {
      j.at("class").get_to(c.device_class);

      j.at("id").get_to(c.disambiguation_string);

      std::stringstream error_context;
      error_context << "\n\tSee the \"" << c.device_class
                    << "\" device with \"id\": \"" << c.disambiguation_string
                    << "\" in the \"devices\" list of your JSON file.";

      for (auto& raw_parameter : j.at("params").items())
      {
        auto key = magic_enum::enum_cast<Parameters>(raw_parameter.key());
        if (key.has_value())
        {
          // NOTE: this is necessary because it doesn't seem like nlohmann/json
          //       handles std::variant out of the box
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
            Log::error() << "\n\tInvalid initial parameter value type: "
                         << "\"" << raw_parameter.key() << "\": "
                         << raw_parameter.value()
                         << " (typed as \"" << raw_parameter.value().type_name()
                         << "\")." << error_context.str() << std::endl;
          }
        }
        else
        {
          Log::error() << "\n\tInitial parameter \"" << raw_parameter.key()
                       << "\" has no value." << error_context.str()
                       << std::endl;
        }
      }

      for (auto& raw_port : j.at("ports").items())
      {
        auto key = magic_enum::enum_cast<Ports>(raw_port.key());
        if (key.has_value())
        {
          const auto& port_value = raw_port.value();
          auto        fail       = [&](const std::string& message)
          {
            std::stringstream ss;
            ss << "\n\t" << message << error_context.str();
            Log::error() << ss.str() << std::endl;
            throw std::runtime_error(ss.str());
          };

          if (port_value.is_object())
          {
            if (!port_value.contains("signal"))
            {
              fail("Signal port object must define \"signal\".");
            }
            port_value.at("signal").get_to(c.ports[key.value()]);

            if (port_value.contains("delay"))
            {
              if (!port_value.at("delay").is_number())
              {
                fail("Signal port \"delay\" must be numeric.");
              }
              const auto delay = port_value.at("delay").template get<RealT>();
              if (delay < 0.0)
              {
                fail("Signal port \"delay\" must be non-negative.");
              }
              if (delay > 0.0)
              {
                c.port_delays[key.value()] = delay;
              }
            }
          }
          else
          {
            port_value.get_to(c.ports[key.value()]);
          }
        }
        else
        {
          Log::error() << "\n\tInvalid port mapping: \"" << raw_port.key()
                       << "\" has no value." << error_context.str()
                       << std::endl;
        }
      }

      if (j.contains("mon"))
      {
        for (auto& raw_monitored_variable : j.at("mon"))
        {
          auto var_name  = raw_monitored_variable.get<std::string>();
          auto monitored = magic_enum::enum_cast<MonitorableVariables>(var_name);
          if (monitored.has_value())
          {
            c.monitored_variables.insert(monitored.value());
          }
          else
          {
            Log::error() << "\n\tInvalid monitored variable: \"" << var_name
                         << "\" in \"mon\" list." << error_context.str()
                         << std::endl;
          }
        }
      }
    }
  } // namespace PhasorDynamics
} // namespace GridKit
