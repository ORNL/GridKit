#pragma once

#include <sstream>
#include <stdexcept>

#include <nlohmann/json.hpp>

#include <GridKit/Model/EMT/Component/Bus/BusDataJSONParser.hpp>
#include <GridKit/Model/EMT/Component/Line/LineLumped/LineLumpedDataJSONParser.hpp>
#include <GridKit/Model/EMT/Component/Load/LoadZ/LoadZDataJSONParser.hpp>
#include <GridKit/Model/EMT/Component/Source/VoltageSource/VoltageSourceDataJSONParser.hpp>
#include <GridKit/Model/EMT/ComponentDataJSONParser.hpp>
#include <GridKit/Model/EMT/Signal/SignalDataJSONParser.hpp>
#include <GridKit/Model/EMT/SystemModelData.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace EMT
  {
    using json          = nlohmann::json;
    using Log           = ::GridKit::Utilities::Logger;
    using MonitorFormat = ::GridKit::Model::VariableMonitorFormat;

    /// JSON parser function implementation for the `SystemModelData` type
    ///
    /// See the `INPUT_FORMAT.md` in `GridKit/Model/EMT` for more information
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

      /// Gets all signals (allows for systems without signals)
      if (j.contains("signals"))
      {
        j.at("signals").get_to(sm.signal);
      }

      /// Gets all components
      for (auto& raw_component : j.at("devices"))
      {
        auto kind = raw_component.at("class").get<std::string>();
        if (kind == "VoltageSource")
        {
          typename SystemModelData<RealT, IdxT>::VoltageSourceDataT source;
          raw_component.get_to(source);
          sm.voltage_source.push_back(source);
        }
        else if (kind == "Machine")
        {
          typename SystemModelData<RealT, IdxT>::MachineDataT machine;
          raw_component.get_to(machine);
          sm.machine.push_back(machine);
        }
        else if (kind == "LineLumped")
        {
          typename SystemModelData<RealT, IdxT>::LineLumpedDataT line;
          raw_component.get_to(line);
          sm.line_lumped.push_back(line);
        }
        else if (kind == "LoadZ")
        {
          typename SystemModelData<RealT, IdxT>::LoadZDataT loadz;
          raw_component.get_to(loadz);
          sm.loadz.push_back(loadz);
        }
        else if (kind == "Switch")
        {
          typename SystemModelData<RealT, IdxT>::SwitchDataT sw;
          raw_component.get_to(sw);
          sm.sw.push_back(sw);
        }
        else if (kind == "Tgov1")
        {
          typename SystemModelData<RealT, IdxT>::Tgov1DataT gov;
          raw_component.get_to(gov);
          sm.gov.push_back(gov);
        }
        else
        {
          Log::error() << "\n\tInvalid device class: \"" << kind << "\". "
                       << "\n\tSee the \"devices\" list in your JSON file."
                       << std::endl;
          throw std::runtime_error("JSON parser failed");
        }
      }
    }
  } // namespace EMT
} // namespace GridKit
