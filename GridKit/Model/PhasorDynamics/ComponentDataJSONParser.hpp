#pragma once

#include <sstream>
#include <string>
#include <type_traits>

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
              typename Terminals,
              typename InputPorts,
              typename OutputPorts,
              typename MonitorableVariables>
      requires std::is_enum_v<Parameters>
               && std::is_enum_v<Terminals>
               && std::is_enum_v<InputPorts>
               && std::is_enum_v<OutputPorts>
               && std::is_enum_v<MonitorableVariables>
    void from_json(const json&                          j,
                   ComponentData<RealT,
                                 IdxT,
                                 Parameters,
                                 Terminals,
                                 InputPorts,
                                 OutputPorts,
                                 MonitorableVariables>& c)
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

      // WARNING --- TEMPORARY --- : the C++ model data is
      // intentionally split into terminals, input ports, and output ports.
      // This backward-compatible flat JSON "ports" parsing SHOULD NOT become
      // the supported design; it only avoids changing every tracked case in
      // this PR. The JSON schema should be split once the ingress changes are merged
      for (auto& raw_port : j.at("ports").items())
      {
        auto terminal = magic_enum::enum_cast<Terminals>(raw_port.key());
        if (terminal.has_value() && terminal.value() != Terminals::SIZE)
        {
          raw_port.value().get_to(c.terminals[terminal.value()]);
          continue;
        }

        auto input_port = magic_enum::enum_cast<InputPorts>(raw_port.key());
        if (input_port.has_value() && input_port.value() != InputPorts::SIZE)
        {
          raw_port.value().get_to(c.input_ports[input_port.value()]);
          continue;
        }

        auto output_port = magic_enum::enum_cast<OutputPorts>(raw_port.key());
        if (output_port.has_value() && output_port.value() != OutputPorts::SIZE)
        {
          raw_port.value().get_to(c.output_ports[output_port.value()]);
          continue;
        }

        Log::error() << "\n\tInvalid port mapping: \"" << raw_port.key()
                     << "\" has no value." << error_context.str()
                     << std::endl;
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
