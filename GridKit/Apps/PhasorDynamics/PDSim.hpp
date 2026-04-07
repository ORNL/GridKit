#pragma once

#include <filesystem>
#include <fstream>
#include <map>
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
      fs::path                           system_model_file;
      double                             dt;
      double                             tmax;
      std::vector<SystemEvent>           events;
      std::vector<Cue>                   schedule;
      fs::path                           output_file;
      fs::path                           reference_file;
      double                             error_tol;
      SystemModelData<>                  model_data;
      std::map<std::string, SystemEvent> event_map; // event id → event
      std::map<std::string, size_t>      fault_map; // event id → fault index
    };

    void from_json(const json& j, StudyData& c)
    {
      j.at("system_model_file").get_to(c.system_model_file);
      j.at("dt").get_to(c.dt);
      j.at("tmax").get_to(c.tmax);

      for (auto& raw_event : j.at("events"))
      {
        auto& ev = c.events.emplace_back();
        raw_event.at("id").get_to(ev.id);
        raw_event.at("type").get_to(ev.type);
        ev.params = raw_event;
      }

      for (auto& raw_action : j.at("schedule"))
      {
        auto& sa = c.schedule.emplace_back();
        raw_action.at("time").get_to(sa.time);
        raw_action.at("event").get_to(sa.event);
        raw_action.at("action").get_to(sa.action);
      }

      if (j.contains("output_file"))
      {
        j.at("output_file").get_to(c.output_file);
        if (!c.output_file.is_absolute())
        {
          // c.output_file =
        }
      }

      if (j.contains("reference_file"))
      {
        j.at("reference_file").get_to(c.reference_file);
      }

      c.error_tol = j.value("error_tolerance", 1.0e-4);
    }

    std::ifstream openFile(const fs::path& file_path)
    {
      if (!exists(file_path))
      {
        Log::error() << "File not found: " << file_path << std::endl;
      }
      auto fs = std::ifstream(file_path);
      if (!fs)
      {
        Log::error() << "Failed to open file: " << file_path << std::endl;
      }
      return fs;
    }

    StudyData parseStudyData(const fs::path& file_path)
    {
      auto data = StudyData(json::parse(openFile(file_path)));

      auto loc = file_path.parent_path();
      if (!data.system_model_file.is_absolute())
      {
        data.system_model_file = loc / data.system_model_file;
      }
      if (!data.reference_file.empty())
      {
        if (!data.reference_file.is_absolute())
        {
          data.reference_file = loc / data.reference_file;
        }
      }

      auto csv        = ::GridKit::Model::VariableMonitorFormat::CSV;
      data.model_data = parseSystemModelData(data.system_model_file);

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
      std::string model_output_file;
      for (const auto& sink : data.model_data.monitor_sink)
      {
        if (sink.format == csv && sink.delim == ",")
        {
          model_output_file = sink.file_name;
        }
      }
      if (model_output_file.empty())
      {
        data.model_data.monitor_sink.emplace_back(data.output_file, csv);
      }
      else if (data.output_file != model_output_file)
      {
        if (exists(data.output_file))
        {
          if ((!is_symlink(data.output_file)) || (read_symlink(data.output_file) != model_output_file))
          {
            Log::error() << "Study output file not usable" << std::endl;
          }
        }
        else
        {
          fs::create_symlink(model_output_file, data.output_file);
        }
      }

      return data;
    }
  } // namespace PhasorDynamics
} // namespace GridKit
