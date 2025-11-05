#pragma once

#include <sstream>
#include <stdexcept>

#include <nlohmann/json.hpp>

#include <GridKit/Model/PhasorDynamics/Bus/BusData.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    using json = nlohmann::json;

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
            std::stringstream msg;
            msg << "\n\tInvalid initial parameter \"" << raw_parameter.key()
                << "\" in \"init\" section." << error_context.str();
            throw std::runtime_error(msg.str());
          }
        }
      }

      j.at("number").get_to(bd.bus_id);

      auto string_class = j.at("class").get<std::string>();
      if (string_class == "bus")
      {
        bd.bus_type = BusData<RealT, IdxT>::BusType::DEFAULT;
      }
      else if (string_class == "infinite_bus")
      {
        bd.bus_type = BusData<RealT, IdxT>::BusType::SLACK;
      }
      else
      {
        std::stringstream msg;
        msg << "\n\tInvalid bus class: \"" << string_class << "\"."
            << error_context.str();
        throw std::runtime_error(msg.str());
      }

      j.at("v_base").get_to(bd.v_base);

      if (j.contains("freq_base"))
      {
        j.at("freq_base").get_to(bd.freq_base);
      }

      if (j.contains("va_base"))
      {
        j.at("va_base").get_to(bd.va_base);
      }

      if (j.contains("mon"))
      {
        for (auto& raw_monitored_variable : j.at("mon"))
        {
          auto monitored = raw_monitored_variable.get<std::string>();
          if (monitored == "Vr")
          {
            bd.monitored_variables.set(static_cast<size_t>(
                BusData<RealT, IdxT>::MonitorableVariables::VR));
          }
          else if (monitored == "Vi")
          {
            bd.monitored_variables.set(static_cast<size_t>(
                BusData<RealT, IdxT>::MonitorableVariables::VI));
          }
          else if (monitored == "Vm")
          {
            bd.monitored_variables.set(static_cast<size_t>(
                BusData<RealT, IdxT>::MonitorableVariables::VM));
          }
          else if (monitored == "Va")
          {
            bd.monitored_variables.set(static_cast<size_t>(
                BusData<RealT, IdxT>::MonitorableVariables::VA));
          }
          else
          {
            std::stringstream msg;
            msg << "\n\tInvalid monitored variable: \"" << monitored
                << "\" in \"mon\" list." << error_context.str();
            throw std::runtime_error(msg.str());
          }
        }
      }
    }
  } // namespace PhasorDynamics
} // namespace GridKit
