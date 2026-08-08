#pragma once

#include <cstddef>
#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <magic_enum/magic_enum.hpp>
#include <nlohmann/json.hpp>

#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModelDataJSONParser.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace PhasorDynamics
  {
    namespace fs = ::std::filesystem;

    using Log = GridKit::Utilities::Logger;

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
     * @brief One monitor output file, what it contains, and what validates it
     *
     * Pairing the sink with its reference keeps a study's outputs and its
     * checks in one place; a spec with no reference file is written but not
     * validated.
     */
    struct MonitorSpec
    {
      /// Path to output file
      fs::path                                output_file;
      /// Output format
      ::GridKit::Model::VariableMonitorFormat format{
          ::GridKit::Model::VariableMonitorFormat::CSV};
      /// Delimiter (CSV only)
      std::string              delim{","};
      /// Column name globs; empty selects every monitored column
      std::vector<std::string> include;
      /// Path to reference file for validation, or empty for none
      fs::path                 reference_file;
      /// Error tolerance (between output file and reference file)
      std::vector<double>      error_tol;
      /// Type of total error (relative or absolute)
      Testing::ErrorType       error_type{Testing::ErrorType::RELATIVE};
      /// Smallest value at which to scale for relative error
      double                   abs_err_threshold{
          Testing::DEFAULT_ABS_ERROR_THRESHOLD};
    };

    /**
     * @brief Data defined in JSON file for parameterized study
     */
    struct StudyData
    {
      /// path to system model JSON file
      fs::path                 system_model_file;
      /// monitor output time step size, or 0 for no intermediate monitoring
      double                   dt_monitor;
      /// max time
      double                   tmax;
      /// relative tolerance for the solver
      double                   rel_tol;
      /// absolute tolerance for the solver
      double                   abs_tol;
      /// fixed solver time step size, or 0 for adaptive stepping
      double                   dt_fixed;
      /// maximum number of solver time steps, or 0 for the IDA default
      std::size_t              max_steps;
      /// set of system events
      std::vector<SystemEvent> events;
      /// path to output file
      fs::path                 output_file;
      /// path to reference file for validation
      fs::path                 reference_file;
      /// Error tolerance (between output file and reference file)
      std::vector<double>      error_tol;
      /// Type of total error (relative or absolute)
      Testing::ErrorType       error_type;
      /// Smallest value at which to scale for relative error
      double                   abs_err_threshold;
      /// Monitor output files declared by the study
      std::vector<MonitorSpec> monitors;
      /// Instance of model data
      SystemModelData<>        model_data;
    };

    using json = ::nlohmann::json;
    using Log  = ::GridKit::Utilities::Logger;

    inline constexpr double DEFAULT_SOLVER_REL_TOL   = 1.0e-7;
    inline constexpr double DEFAULT_SOLVER_ABS_TOL   = 1.0e-9;
    inline constexpr double DEFAULT_VERIFICATION_TOL = 1.0e-4;

    /**
     * @brief Read the error settings shared by a study and its monitors
     */
    inline void parseErrorSettings(const json&          j,
                                   std::vector<double>& error_tol,
                                   Testing::ErrorType&  error_type,
                                   double&              abs_err_threshold)
    {
      using ErrorType = Testing::ErrorType;

      if (j.contains("error_tolerance"))
      {
        const auto& tolj = j.at("error_tolerance");
        if (tolj.is_array())
        {
          tolj.get_to(error_tol);
        }
        else
        {
          tolj.get_to(error_tol.emplace_back());
        }
      }
      else
      {
        error_tol.push_back(DEFAULT_VERIFICATION_TOL);
      }

      if (j.contains("error_type"))
      {
        auto type_str  = j.at("error_type").get<std::string>();
        auto type_wrap = magic_enum::enum_cast<ErrorType>(
            type_str, magic_enum::case_insensitive);
        if (!type_wrap.has_value())
        {
          Log::error() << "Invalid error type \"" << type_str << "\"; "
                       << "must be either \"relative\" or \"absolute\"";
        }
        error_type = type_wrap.value();
      }
      else
      {
        error_type = ErrorType::RELATIVE;
      }

      abs_err_threshold =
          j.value("abs_err_threshold", Testing::DEFAULT_ABS_ERROR_THRESHOLD);
    }

    /**
     * @brief JSON parser implemntation for `StudyData`
     */
    void from_json(const json& j, StudyData& c)
    {
      using namespace magic_enum;

      j.at("system_model_file").get_to(c.system_model_file);
      c.dt_monitor = j.value("dt_monitor", 0.0);
      j.at("tmax").get_to(c.tmax);
      c.rel_tol   = j.value("rel_tol", DEFAULT_SOLVER_REL_TOL);
      c.abs_tol   = j.value("abs_tol", DEFAULT_SOLVER_ABS_TOL);
      c.dt_fixed  = j.value("dt_fixed", 0.0);
      c.max_steps = j.value("max_steps", std::size_t{0});

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

      parseErrorSettings(j, c.error_tol, c.error_type, c.abs_err_threshold);

      if (j.contains("monitors"))
      {
        for (const auto& raw_mon : j.at("monitors"))
        {
          auto& mon = c.monitors.emplace_back();
          raw_mon.at("file_name").get_to(mon.output_file);
          mon.delim   = raw_mon.value("delim", std::string{","});
          mon.include = parseIncludeList(raw_mon);

          if (raw_mon.contains("format"))
          {
            auto fmt_str = raw_mon.at("format").get<std::string>();
            auto fmt     = enum_cast<::GridKit::Model::VariableMonitorFormat>(
                fmt_str, case_insensitive);
            if (!fmt.has_value())
            {
              Log::error() << "Invalid monitor format \"" << fmt_str << "\"."
                           << std::endl;
            }
            mon.format = fmt.value();
          }

          if (raw_mon.contains("reference_file"))
          {
            raw_mon.at("reference_file").get_to(mon.reference_file);
          }

          parseErrorSettings(raw_mon, mon.error_tol, mon.error_type, mon.abs_err_threshold);
        }
      }
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

      // Monitors declared by the study become sinks on the model. Reference
      // paths are study-relative like `reference_file`; output paths are not,
      // so they land in the working directory as they always have.
      for (auto& mon : data.monitors)
      {
        if (!mon.reference_file.empty() && !mon.reference_file.is_absolute())
        {
          mon.reference_file = loc / mon.reference_file;
        }
        data.model_data.monitor_sink.emplace_back(
            mon.format, mon.output_file.string(), mon.delim, mon.include);
      }

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
        if (!data.output_file.empty())
        {
          // Add study output file to model if one did not already exist
          data.model_data.monitor_sink.emplace_back(csv, data.output_file);
        }
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

    void checkCommandLine(int argc, const std::string& appName)
    {
      if (argc < 2)
      {
        Log::error() << "No input file provided" << std::endl;
        std::cout << std::format(
            "\n"
            "Usage:\n"
            "       {} <json-input-file>\n"
            "\n"
            "Please provide a json input file for the study to run.\n"
            "\n",
            appName);
        exit(1);
      }
    }

    Testing::TestStatus checkErrors(
        const StudyData& study_data,
        bool             print_results = true)
    {
      // Generate aggregate errors comparing variable output to reference solution
      auto func   = std::string{"monitor file vs reference file"};
      auto status = Testing::TestStatus{func.c_str()};

      auto compare = [&status, print_results](const fs::path&            out_file,
                                              const fs::path&            ref_file,
                                              Testing::ErrorType         error_type,
                                              double                     abs_err_threshold,
                                              const std::vector<double>& error_tol)
      {
        if (out_file.empty() || ref_file.empty())
        {
          return;
        }

        auto errorSet = Testing::compareCSV(out_file,
                                            ref_file,
                                            error_type,
                                            abs_err_threshold);

        // Print the errors
        if (print_results)
        {
          std::cout << out_file.string() << " vs " << ref_file.string() << '\n';
          errorSet->display();
        }

        // Check against specified tolerance
        if (error_tol.size() > 1)
        {
          status *= error_tol.size() == errorSet->var_errors.size();
          for (std::size_t i = 0; i < error_tol.size(); ++i)
          {
            status *= errorSet->var_errors[i].max_value < error_tol[i];
          }
        }
        else
        {
          status *= errorSet->total_error.max_value < error_tol[0];
        }
      };

      compare(study_data.output_file,
              study_data.reference_file,
              study_data.error_type,
              study_data.abs_err_threshold,
              study_data.error_tol);

      for (const auto& mon : study_data.monitors)
      {
        compare(mon.output_file,
                mon.reference_file,
                mon.error_type,
                mon.abs_err_threshold,
                mon.error_tol);
      }

      if (print_results)
      {
        status.report();
      }
      return status;
    }

  } // namespace PhasorDynamics
} // namespace GridKit
