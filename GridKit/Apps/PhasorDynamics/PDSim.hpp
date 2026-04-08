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

    /**
     * @brief Describes an event that is used to modify the simulation at the
     * given time point
     */
    struct SystemEvent
    {
      /// Type of event determines action performed
      enum class Type
      {
        FAULT_ON,
        FAULT_OFF
      };

      /// Time event takes place
      double      time;
      /// Event type
      Type        type;
      /// ID of element used in event (e.g., bus fault id)
      std::size_t element_id;
    };

    /**
     * @brief Data defined in JSON file for parameterized study
     */
    struct StudyData
    {
      /// path to system model JSON file
      fs::path                 system_model_file;
      /// time step size
      double                   dt;
      /// max time
      double                   tmax;
      /// set of system events
      std::vector<SystemEvent> events;
      /// path to output file
      fs::path                 output_file;
      /// path to reference file for validation
      fs::path                 reference_file;
      /// Error tolerance (between output file and reference file)
      double                   error_tol;
      /// Instance of model data
      SystemModelData<>        model_data;
    };

    using json = ::nlohmann::json;
    using Log  = ::GridKit::Utilities::Logger;

    /**
     * @brief JSON parser implemntation for `StudyData`
     */
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
      }

      if (j.contains("reference_file"))
      {
        j.at("reference_file").get_to(c.reference_file);
      }

      c.error_tol = j.value("error_tolerance", 1.0e-4);
    }

    /**
     * @brief Check for existence and successful input file open
     */
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

    /**
     * @brief Wrapper function to parse `StudyData` from JSON and perform
     * follow-up configuration
     */
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
      // Find output file (CSV) specified in model input file
      for (const auto& sink : data.model_data.monitor_sink)
      {
        if (sink.format == csv && sink.delim == ",")
        {
          model_output_file = sink.file_name;
        }
      }

      if (model_output_file.empty())
      {
        // Add study output file to model if one did not already exist
        data.model_data.monitor_sink.emplace_back(data.output_file, csv);
      }
      else
      {
        if (data.output_file.empty())
        {
          data.output_file = model_output_file;
        }
        else
        {
          // If model file already specifies a CSV output file, then the study
          // output file must be a symlink to the model output file
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
      }

      return data;
    }
  } // namespace PhasorDynamics
} // namespace GridKit
