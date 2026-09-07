#pragma once

#include <stdexcept>

#include <nlohmann/json.hpp>

#include <GridKit/Model/EMT/ContainerDataJSONParser.hpp>
#include <GridKit/Model/EMT/SystemModelData.hpp>

namespace GridKit
{
  namespace EMT
  {
    using json          = nlohmann::json;
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

      validateJsonFields(j, "System model", {"header", "monitors", "signals", "devices", "inputs", "outputs"});
      const auto& header = j.at("header");
      validateJsonFields(header, "System model header", {"format_version", "format_revision", "case_name", "case_date_time", "case_description", "case_comments"});

      if (header.contains("format_version"))
      {
        sm.format_version = parseFiniteReal<double>(header.at("format_version"), "format_version");
      }

      if (header.contains("format_revision"))
      {
        const auto revision = parseFiniteReal<double>(header.at("format_revision"), "format_revision");
        if (revision < 0 || revision > std::numeric_limits<unsigned short>::max() || std::trunc(revision) != revision)
          throw std::invalid_argument("format_revision requires a nonnegative integer within the revision range");
        sm.format_revision = static_cast<unsigned short>(revision);
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
        if (!j.at("monitors").is_array())
          throw std::invalid_argument("System model monitors must be an array");
        for (const auto& raw_mon : j.at("monitors"))
        {
          validateJsonFields(raw_mon, "Monitor sink", {"file_name", "format", "delim"});
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
            throw std::invalid_argument("Invalid monitor format \"" + fmt_str + "\"");
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
