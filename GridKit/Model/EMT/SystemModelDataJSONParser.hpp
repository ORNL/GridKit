#pragma once

#include <stdexcept>

#include <nlohmann/json.hpp>

#include <GridKit/Model/EMT/ContainerDataJSONParser.hpp>
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

      parseContainerData(j, static_cast<ContainerData<RealT, IdxT>&>(sm), "root");
      if (!sm.inputs.empty())
      {
        throw std::runtime_error("The root SystemModel cannot bind Container inputs");
      }
    }
  } // namespace EMT
} // namespace GridKit
