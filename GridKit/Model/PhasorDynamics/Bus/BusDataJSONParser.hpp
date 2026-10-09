#pragma once

#include <sstream>
#include <stdexcept>

#include <magic_enum/magic_enum.hpp>
#include <nlohmann/json.hpp>

#include <GridKit/Model/PhasorDynamics/Bus/BusData.hpp>
#include <GridKit/Utilities/Enum.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    using json = nlohmann::json;
    using Log  = ::GridKit::Utilities::Logger;

    /// JSON parser function implementation for the `BusData` type
    ///
    /// See the `README.md` in `GridKit/Model/PhasorDynamics` for more information
    template <typename RealT, typename IdxT>
    void from_json(const json& j, BusData<RealT, IdxT>& bd)
    {
      j.at("name").get_to(bd.name);

      std::stringstream error_context;
      error_context << "\n\tSee bus number " << bd.bus_id
                    << " (\"name\": \"" << bd.name << "\") "
                    << "in the \"buses\" list of your JSON file.";

      if (j.contains("init"))
      {
        for (auto& raw_parameter : j.at("init").items())
        {
          if (raw_parameter.key() == "Vi")
          {
            raw_parameter.value().get_to(bd.Vi0);
          }
          else if (raw_parameter.key() == "Vr")
          {
            raw_parameter.value().get_to(bd.Vr0);
          }
          else
          {
            Log::error() << "\n\tInvalid initial parameter \""
                         << raw_parameter.key() << "\" in \"init\" section."
                         << error_context.str() << std::endl;
          }
        }
      }

      j.at("number").get_to(bd.bus_id);

      auto string_class = j.at("class").get<std::string>();
      if (string_class == "Bus")
      {
        bd.bus_type = BusData<RealT, IdxT>::BusType::DEFAULT;
      }
      else if (string_class == "BusInfinite")
      {
        bd.bus_type = BusData<RealT, IdxT>::BusType::SLACK;
      }
      else if (string_class == "SignalVoltageIn")
      {
        bd.bus_type = BusData<RealT, IdxT>::BusType::SIGNAL_VOLTAGE_IN;
      }
      else if (string_class == "SignalVoltageOut")
      {
        bd.bus_type = BusData<RealT, IdxT>::BusType::SIGNAL_VOLTAGE_OUT;
      }
      else
      {
        Log::error() << "\n\tInvalid bus class: \"" << string_class << "\"."
                     << error_context.str() << std::endl;
        throw std::runtime_error("JSON parser failed");
      }

      using Parameters = typename BusData<RealT, IdxT>::Parameters;
      for (auto& raw_parameter : j.at("params").items())
      {
        auto key = Utilities::enum_cast<Parameters>(raw_parameter.key());
        if (key.has_value())
        {
          // NOTE: this is necessary because it doesn't seem like nlohmann/json
          //       handles std::variant out of the box
          if (raw_parameter.value().is_boolean())
          {
            bd.parameters[key.value()] = raw_parameter.value().template get<bool>();
          }
          else if (raw_parameter.value().is_number_float())
          {
            bd.parameters[key.value()] = raw_parameter.value().template get<RealT>();
          }
          else if (raw_parameter.value().is_number_integer())
          {
            bd.parameters[key.value()] = raw_parameter.value().template get<IdxT>();
          }
          else
          {
            Log::error() << "\n\tInvalid bus parameter value type: "
                         << "\"" << raw_parameter.key() << "\": "
                         << raw_parameter.value()
                         << " (typed as \"" << raw_parameter.value().type_name()
                         << "\")." << error_context.str() << std::endl;
            throw std::runtime_error("JSON parser failed");
          }
        }
        else
        {
          Log::error() << "\n\tBus parameter \"" << raw_parameter.key()
                       << "\" has no value." << error_context.str()
                       << std::endl;
        }
      }

      using SignalInputs = typename BusData<RealT, IdxT>::SignalInputs;
      if (j.contains("ports_in"))
      {
        for (auto& raw_input : j.at("ports_in").items())
        {
          auto input = Utilities::enum_cast<SignalInputs>(raw_input.key());
          if (input.has_value()
              && static_cast<size_t>(input.value())
                     < Utilities::enum_size<SignalInputs>())
          {
            raw_input.value().get_to(bd.signal_inputs[input.value()]);
          }
        }
      }

      using SignalOutputs = typename BusData<RealT, IdxT>::SignalOutputs;
      if (j.contains("ports_out"))
      {
        for (auto& raw_output : j.at("ports_out").items())
        {
          auto output = Utilities::enum_cast<SignalOutputs>(raw_output.key());
          if (output.has_value()
              && static_cast<size_t>(output.value())
                     < Utilities::enum_size<SignalOutputs>())
          {
            raw_output.value().get_to(bd.signal_outputs[output.value()]);
          }
        }
      }

      if (j.contains("mon"))
      {
        using MonitorableVariables = typename BusData<RealT, IdxT>::MonitorableVariables;
        for (auto& raw_monitored_variable : j.at("mon"))
        {
          auto var_name  = raw_monitored_variable.get<std::string>();
          auto monitored = Utilities::enum_cast<MonitorableVariables>(var_name, magic_enum::case_insensitive);
          if (monitored.has_value())
          {
            bd.monitored_variables.insert(monitored.value());
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
