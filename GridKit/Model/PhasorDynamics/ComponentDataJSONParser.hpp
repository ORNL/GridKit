#pragma once

#include <sstream>
#include <stdexcept>

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
              typename MonitorableVariables,
              typename ExternalVariables>
      requires std::is_enum_v<Parameters>
               && std::is_enum_v<Ports>
               && std::is_enum_v<MonitorableVariables>
               && std::is_enum_v<ExternalVariables>
               && requires { ExternalVariables::MAXIMUM; }
    void from_json(const json& j, ComponentData<RealT, IdxT, Parameters, Ports, MonitorableVariables, ExternalVariables>& c)
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

      if (j.contains("init"))
      {
        for (auto& raw_initial_value : j.at("init").items())
        {
          using magic_enum::case_insensitive;

          auto key = magic_enum::enum_cast<ExternalVariables>(
              raw_initial_value.key(), case_insensitive);
          if (!key.has_value() || key.value() == ExternalVariables::MAXIMUM)
          {
            Log::error() << "\n\tInvalid external initial value: \""
                         << raw_initial_value.key()
                         << "\" has no value." << error_context.str()
                         << std::endl;
            ++c.input_error_count;
            continue;
          }

          if (c.external_initial_values.contains(key.value()))
          {
            Log::error() << "\n\tDuplicate external initial value: \""
                         << raw_initial_value.key() << "\"."
                         << error_context.str() << std::endl;
            ++c.input_error_count;
            continue;
          }

          const auto& value = raw_initial_value.value();
          if (value.is_boolean())
          {
            c.external_initial_values[key.value()] = value.template get<bool>();
          }
          else if (value.is_number_float())
          {
            c.external_initial_values[key.value()] = value.template get<RealT>();
          }
          else if (value.is_number_integer())
          {
            c.external_initial_values[key.value()] = value.template get<IdxT>();
          }
          else
          {
            Log::error() << "\n\tInvalid initial value type: "
                         << "\"" << raw_initial_value.key() << "\": "
                         << value
                         << " (typed as \"" << value.type_name()
                         << "\")." << error_context.str() << std::endl;
            ++c.input_error_count;
          }
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
