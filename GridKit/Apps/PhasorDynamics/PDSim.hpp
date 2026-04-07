#pragma once

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <map>
#include <stdexcept>
#include <vector>

#include <nlohmann/json.hpp>

#include <GridKit/Model/PhasorDynamics/BusFault/BusFaultData.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace fs = ::std::filesystem;

    using json = ::nlohmann::json;
    using Log  = ::GridKit::Utilities::Logger;

    struct SystemEvent
    {
      std::string id;
      std::string type;   // "bus_fault", future: "branch_trip", etc.
      json        params; // type-specific data
    };

    struct Cue
    {
      double      time;
      std::string event;  // references SystemEvent::id
      std::string action; // "on", "off", etc.
    };

    struct StudyData
    {
      using MonitorSinkSpec = SystemModelData<>::MonitorSinkSpec;

      int                                format_version{1};
      std::string                        name;
      std::string                        description;
      fs::path                           case_file;
      double                             dt;
      double                             tmax;
      std::vector<SystemEvent>           events;
      std::vector<Cue>                   schedule;
      MonitorSinkSpec                    output;
      fs::path                           reference_file;
      double                             error_tol;
      SystemModelData<>                  model_data;
      std::map<std::string, SystemEvent> event_map; // event id → event
      std::map<std::string, size_t>      fault_map; // event id → fault index
    };

    inline void from_json(const json& j, StudyData& c)
    {
      c.format_version = j.value("format_version", 1);
      c.name           = j.value("study_name", std::string{});
      c.description    = j.value("study_description", std::string{});

      j.at("case_file").get_to(c.case_file);
      j.at("dt").get_to(c.dt);
      j.at("tmax").get_to(c.tmax);

      if (j.contains("events"))
      {
        for (auto& raw_event : j.at("events"))
        {
          auto& ev = c.events.emplace_back();
          raw_event.at("id").get_to(ev.id);
          raw_event.at("type").get_to(ev.type);
          ev.params = raw_event.at("params");
        }
      }

      if (j.contains("schedule"))
      {
        for (auto& raw_action : j.at("schedule"))
        {
          auto& sa = c.schedule.emplace_back();
          raw_action.at("time").get_to(sa.time);
          raw_action.at("event").get_to(sa.event);
          raw_action.at("action").get_to(sa.action);
        }
      }

      {
        using Format = ::GridKit::Model::VariableMonitorFormat;
        if (j.contains("output"))
        {
          const auto& out    = j.at("output");
          c.output.file_name = out.value("file_name", std::string{});
          c.output.delim     = out.value("delim", std::string(","));

          auto fmt_str = out.value("format", std::string("CSV"));
          std::transform(fmt_str.begin(), fmt_str.end(), fmt_str.begin(), ::toupper);

          if (fmt_str == "CSV")
            c.output.format = Format::CSV;
          else if (fmt_str == "JSON")
            c.output.format = Format::JSON;
          else if (fmt_str == "YAML")
            c.output.format = Format::YAML;
          else
            throw std::runtime_error("Invalid output format: \"" + fmt_str + "\"");
        }
        else
        {
          c.output.format = Format::CSV;
        }
      }

      if (j.contains("reference_file"))
      {
        j.at("reference_file").get_to(c.reference_file);
      }

      c.error_tol = j.value("error_tol", 1.0e-4);
    }

    inline std::ifstream openFile(const fs::path& file_path)
    {
      if (!exists(file_path))
      {
        Log::error() << "File not found: " << file_path << std::endl;
        throw std::runtime_error("File not found: " + file_path.string());
      }
      auto fs = std::ifstream(file_path);
      if (!fs)
      {
        Log::error() << "Failed to open file: " << file_path << std::endl;
        throw std::runtime_error("Failed to open file: " + file_path.string());
      }
      return fs;
    }

    inline StudyData parseStudyData(const fs::path& file_path)
    {
      auto data = StudyData(json::parse(openFile(file_path)));

      if (data.format_version != 1)
      {
        throw std::runtime_error(
            "Unsupported study file format_version: " + std::to_string(data.format_version));
      }

      auto loc = file_path.parent_path();
      if (!data.case_file.is_absolute())
      {
        data.case_file = loc / data.case_file;
      }
      if (!data.reference_file.empty())
      {
        if (!data.reference_file.is_absolute())
        {
          data.reference_file = loc / data.reference_file;
        }
      }
      if (!data.output.file_name.empty())
      {
        fs::path p(data.output.file_name);
        if (!p.is_absolute())
        {
          data.output.file_name = (loc / p).string();
        }
      }

      data.model_data = parseSystemModelData(data.case_file);

      // Build event lookups and inject event data into model
      size_t fault_idx = 0;
      for (const auto& ev : data.events)
      {
        data.event_map[ev.id] = ev;

        if (ev.type == "bus_fault")
        {
          using BFD = BusFaultData<double, size_t>;
          BFD bfd;
          bfd.parameters[BFD::Parameters::R]      = ev.params.value("R", 0.0);
          bfd.parameters[BFD::Parameters::X]      = ev.params.value("X", 1e-5);
          bfd.parameters[BFD::Parameters::state0] = false;
          bfd.ports[BFD::Ports::bus]              = ev.params.at("bus").get<size_t>();
          data.model_data.bus_fault.push_back(bfd);
          data.fault_map[ev.id] = fault_idx++;
        }
      }

      // Validate schedule
      for (size_t i = 0; i < data.schedule.size(); ++i)
      {
        const auto& cue = data.schedule[i];

        if (data.event_map.find(cue.event) == data.event_map.end())
        {
          throw std::runtime_error(
              "Schedule cue at time " + std::to_string(cue.time)
              + " references undefined event: \"" + cue.event + "\"");
        }

        if (cue.action != "on" && cue.action != "off")
        {
          throw std::runtime_error(
              "Schedule cue at time " + std::to_string(cue.time)
              + " has unrecognized action: \"" + cue.action + "\"");
        }

        if (i > 0 && cue.time < data.schedule[i - 1].time)
        {
          throw std::runtime_error(
              "Schedule is not sorted by time: cue at "
              + std::to_string(cue.time) + " follows "
              + std::to_string(data.schedule[i - 1].time));
        }
      }

      data.model_data.monitor_sink.clear();
      data.model_data.monitor_sink.push_back(data.output);

      return data;
    }
  } // namespace PhasorDynamics
} // namespace GridKit
