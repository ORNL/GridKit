#pragma once

#include <cassert>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include <GridKit/Model/PhasorDynamics/Component.hpp>
#include <GridKit/Model/PhasorDynamics/ComponentDataJSONParser.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace fs = ::std::filesystem;

    using json = ::nlohmann::json;
    using Log  = ::GridKit::Utilities::Logger;

    /// One scheduled cue: at `time`, deliver `action` to `target`.
    /// Consumed by PDSim — `time` drives integration, `target`+`action` are
    /// forwarded to `SystemModel::cue`. Components themselves never see this.
    struct Cue
    {
      double      time{0.0};
      std::string target;
      Action      action{Action::Off};
    };

    struct StudyData
    {
      int               format_version{1};
      fs::path          case_file;
      double            dt{0.0};
      double            tmax{0.0};
      std::vector<Cue>  schedule;
      fs::path          output_file;
      fs::path          reference_file;
      double            error_tol{1.0e-4};
      SystemModelData<> model_data;
    };

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
      auto raw = json::parse(openFile(file_path));

      StudyData data;

      data.format_version = raw.value("format_version", 1);
      if (data.format_version != 1)
      {
        throw std::runtime_error(
            "Unsupported solver file format_version: "
            + std::to_string(data.format_version));
      }

      raw.at("case_file").get_to(data.case_file);
      auto loc = file_path.parent_path();
      if (!data.case_file.is_absolute())
      {
        data.case_file = loc / data.case_file;
      }

      raw.at("dt").get_to(data.dt);
      raw.at("tmax").get_to(data.tmax);

      if (raw.contains("output_file"))
      {
        raw.at("output_file").get_to(data.output_file);
      }
      if (raw.contains("reference_file"))
      {
        raw.at("reference_file").get_to(data.reference_file);
        if (!data.reference_file.is_absolute())
        {
          data.reference_file = loc / data.reference_file;
        }
      }
      data.error_tol = raw.value("error_tolerance", 1.0e-4);

      data.model_data = parseSystemModelData(data.case_file);

      if (raw.contains("schedule"))
      {
        double last_time = -std::numeric_limits<double>::infinity();
        for (const auto& raw_cue : raw.at("schedule"))
        {
          Cue cue;
          raw_cue.at("time").get_to(cue.time);
          raw_cue.at("target").get_to(cue.target);

          std::string action_str;
          raw_cue.at("action").get_to(action_str);
          if (action_str == "on")
            cue.action = Action::On;
          else if (action_str == "off")
            cue.action = Action::Off;
          else
            throw std::runtime_error(
                "Solver file: unknown action '" + action_str + "' (valid: on, off)");

          assert(cue.time <= data.tmax);
          if (cue.time < last_time)
          {
            throw std::runtime_error(
                "Solver file: schedule is not monotone non-decreasing in time");
          }
          last_time = cue.time;
          data.schedule.push_back(std::move(cue));
        }
      }

      auto        csv = ::GridKit::Model::VariableMonitorFormat::CSV;
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
        if (!data.output_file.empty())
        {
          data.model_data.monitor_sink.emplace_back(data.output_file.string(), csv);
        }
      }
      else
      {
        if (data.output_file.empty())
        {
          data.output_file = model_output_file;
        }
        else if (data.output_file.string() != model_output_file)
        {
          if (exists(data.output_file))
          {
            if ((!is_symlink(data.output_file))
                || (read_symlink(data.output_file) != model_output_file))
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
