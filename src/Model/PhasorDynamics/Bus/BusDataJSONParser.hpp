#pragma once

#include <Model/PhasorDynamics/Bus/BusData.hpp>
#include <nlohmann/json.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    using json = nlohmann::json;

    /// JSON parser function implementation for the `BusData` type
    ///
    /// See the `README.md` in `src/Model/PhasorDynamics` for more information
    template <typename RealT, typename IdxT>
    void from_json(const json& j, BusData<RealT, IdxT>& bd)
    {
      j.at("name").get_to(bd.name);

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
            throw "Invalid initial parameter";
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
        throw "Invalid bus class";
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
            throw "Invalid monitored variable";
          }
        }
      }
    }
  } // namespace PhasorDynamics
} // namespace GridKit
