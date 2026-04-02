#pragma once

#include <filesystem>
#include <fstream>
#include <vector>

#include <magic_enum/magic_enum.hpp>
#include <nlohmann/json.hpp>

#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace fs = ::std::filesystem;

    struct SystemEvent
    {
      enum class Type
      {
        FAULT_ON,
        FAULT_OFF
      };

      double      time;
      Type        type;
      std::size_t element_id;
    };

    struct StudyData
    {
      fs::path                 system_model_file;
      double                   dt;
      double                   tmax;
      std::vector<SystemEvent> events;
      fs::path                 output_file;
      fs::path                 reference_file;
      double                   error_tol;
      SystemModelData<>        model_data;
    };

    using json = ::nlohmann::json;
    using Log  = ::GridKit::Utilities::Logger;

    void from_json(const json& j, StudyData& c)
    {
      using namespace magic_enum;

      j.at("system_model_file").get_to(c.system_model_file);
      j.at("dt").get_to(c.dt);
      j.at("tmax").get_to(c.tmax);

      for (auto& raw_event : j.at("events"))
      {
        auto& event = c.events.emplace_back();
        raw_event.at("time").get_to(event.time);
        raw_event.at("element_id").get_to(event.element_id);

        auto type_str   = raw_event.at("type").get<std::string>();
        using EventType = SystemEvent::Type;
        auto type_wrap  = enum_cast<EventType>(type_str, case_insensitive);
        if (!type_wrap.has_value())
        {
          Log::error() << "Unable to parse event type \"" << type_str << "\"\n";
        }
        event.type = type_wrap.value();
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
      else
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
