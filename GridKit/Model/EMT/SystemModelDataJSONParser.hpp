/**
 * @file SystemModelDataJSONParser.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief JSON parsing of EMT cases.
 */

#pragma once

#include <string>

#include <nlohmann/json.hpp>

#include <GridKit/Model/EMT/SystemModelData.hpp>

namespace GridKit
{
  namespace EMT
  {
    using json = nlohmann::json;

    /// The caller sets the id
    template <typename IdxT>
    void from_json(const json& j, DeviceData<IdxT>& device)
    {
      j.at("class").get_to(device.device_class);
      device.params = j.value("params", json::object());

      // A number names a bus, and a string a variable or an output as `owner.name`
      const auto references = [&j](const char* section)
      {
        std::map<std::string, std::string> named;
        const json                         entries = j.value(section, json::object());
        for (const auto& [name, value] : entries.items())
        {
          if (value.is_string())
          {
            named[name] = value.template get<std::string>();
          }
          else
          {
            named[name] = std::to_string(value.template get<IdxT>());
          }
        }
        return named;
      };
      device.terminals = references("terminals");
      device.inputs    = references("inputs");
      device.outputs   = references("outputs");
      device.mon       = j.value("mon", device.mon);
    }

    template <typename RealT, typename IdxT>
    void from_json(const json& j, SystemModelData<RealT, IdxT>& data)
    {
      j.at("params").at("freq_base").get_to(data.freq_base);
      j.at("params").at("va_base").get_to(data.va_base);

      for (const auto& entry : j.at("buses"))
      {
        auto& bus = data.buses.emplace_back(entry.get<DeviceData<IdxT>>());
        bus.id    = std::to_string(entry.at("number").get<IdxT>());
      }
      for (const auto& entry : j.at("devices"))
      {
        auto& device = data.devices.emplace_back(entry.get<DeviceData<IdxT>>());
        entry.at("id").get_to(device.id);
      }
    }
  } // namespace EMT
} // namespace GridKit
