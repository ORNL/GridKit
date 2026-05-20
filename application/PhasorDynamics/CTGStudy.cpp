#include "CTGStudy.hpp"

#include <string>

#include <magic_enum/magic_enum.hpp>
#include <nlohmann/json.hpp>

#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    using namespace GridKit::Testing;

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

    CTGStudy::CTGStudy(const StudyData& data)
      : study_data_(data),
        sys_(study_data_.model_data),
        ida_(&sys_)
    {
      sys_.allocate();
      ida_.configureSimulation();
    }

    void CTGStudy::run()
    {
      using EventType = SystemEvent::Type;
      using Clock     = std::chrono::high_resolution_clock;

      // Start timer
      const auto start = Clock::now();

      // Initilize simultation for first run
      RealT dt         = study_data_.dt;
      RealT final_time = study_data_.tmax;
      RealT curr_time  = 0.0;
      ida_.initializeSimulation(0.0, false);
      for (const auto& event : study_data_.events)
      {
        // Run to event time
        int nout = static_cast<int>(std::round((event.time - curr_time) / dt));
        ida_.runSimulation(event.time, nout);

        // Set up run for event (to start at event time)
        switch (event.type)
        {
        case EventType::FAULT_ON:
          sys_.getBusFault(event.element_id)->setStatus(true);
          break;
        case EventType::FAULT_OFF:
          sys_.getBusFault(event.element_id)->setStatus(false);
          break;
        }

        // Re-initialize simulation at event time
        ida_.initializeSimulation(event.time, true);
        curr_time = event.time;
      }

      // Run to final time
      int nout = static_cast<int>(std::round((final_time - curr_time) / dt));
      ida_.runSimulation(final_time, nout);

      // Stop the variable monitor
      sys_.stopMonitor();

      // Save run time
      const auto stop = Clock::now();
      dur_            = stop - start;
    }

    TestStatus CTGStudy::checkStatus(Print print) const
    {
      // Generate aggregate errors comparing variable output to reference solution
      std::string func{"monitor file vs reference file"};
      TestStatus  status{func.c_str()};
      const auto& out_file = study_data_.output_file;
      const auto& ref_file = study_data_.reference_file;
      if (!out_file.empty() && !ref_file.empty())
      {
        auto errorSet = compareCSV(out_file, ref_file);

        // Print the errors
        if (print)
        {
          errorSet.display();
        }

        // Check against specified tolerance
        status *= errorSet.total.max_error < study_data_.error_tol;

        if (print)
        {
          status.report();
        }
      }
      return status;
    }

  } // namespace PhasorDynamics
} // namespace GridKit
