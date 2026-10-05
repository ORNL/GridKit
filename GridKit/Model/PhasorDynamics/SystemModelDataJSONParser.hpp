#pragma once

#include <sstream>
#include <stdexcept>
#include <string_view>

#include <nlohmann/json.hpp>

#include <GridKit/Model/PhasorDynamics/Bus/BusDataJSONParser.hpp>
#include <GridKit/Model/PhasorDynamics/ComponentDataJSONParser.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNodeDataJSONParser.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    using json          = nlohmann::json;
    using Log           = ::GridKit::Utilities::Logger;
    using MonitorFormat = ::GridKit::Model::VariableMonitorFormat;

    /// JSON parser function implementation for the `SystemModelData` type
    ///
    /// See the `README.md` in `GridKit/Model/PhasorDynamics` for more information
    template <typename RealT = double, typename IdxT = size_t>
    void from_json(const json& j, SystemModelData<RealT, IdxT>& sm)
    {
      auto enum_parse = []<typename EnumT, typename KeyT>(EnumT, KeyT&& key)
      {
        return magic_enum::enum_cast<EnumT>(key, magic_enum::case_insensitive);
      };

      auto header = j.at("header");

      if (header.contains("format_version"))
      {
        header.at("format_version").get_to(sm.format_version);
      }

      if (header.contains("format_revision"))
      {
        header.at("format_revision").get_to(sm.format_revision);
      }

      header.at("case_name").get_to(sm.case_name);

      if (header.contains("case_date_time"))
      {
        header.at("case_date_time").get_to(sm.case_date_time);
      }

      header.at("case_description").get_to(sm.case_description);
      header.at("case_comments").get_to(sm.case_comments);

      const auto& params = j.at("params");
      params.at("freq_base").get_to(sm.freq_base);
      params.at("va_base").get_to(sm.va_base);

      if (j.contains("monitors"))
      {
        for (auto&& raw_mon : j.at("monitors"))
        {
          auto file_name = raw_mon.value("file_name", std::string{});
          auto fmt_str   = raw_mon.at("format").get<std::string>();
          auto format    = enum_parse(MonitorFormat{}, fmt_str);
          auto delim     = raw_mon.value("delim", std::string(","));
          if (format.has_value())
          {
            sm.monitor_sink.emplace_back(format.value(), file_name, delim);
          }
          else
          {
            Log::error() << "\n\tInvalid monitor format: \"" << fmt_str << "\"."
                         << "\n\tSee the \"monitors\" list in your JSON file."
                         << std::endl;
          }
        }
      }

      /// Gets all electrical buses
      j.at("buses").get_to(sm.bus);

      /// Gets all signal nodes (allows for systems without signals)
      if (j.contains("signals"))
      {
        j.at("signals").get_to(sm.signal);
      }

      /// Gets all components
      for (auto& raw_component : j.at("devices"))
      {
        const auto kind   = raw_component.at("class").get<std::string>();
        bool       parsed = false;

        forEachDeviceList(sm,
                          [&](std::string_view device_class, auto& devices)
                          {
                            if (kind != device_class)
                            {
                              return;
                            }
                            raw_component.get_to(devices.emplace_back());
                            parsed = true;
                          });

        if (!parsed)
        {
          Log::error() << "\n\tInvalid device class: \"" << kind << "\". "
                       << "\n\tSee the \"devices\" list in your JSON file."
                       << std::endl;
          throw std::runtime_error("JSON parser failed");
        }
      }
    }
  } // namespace PhasorDynamics
} // namespace GridKit
