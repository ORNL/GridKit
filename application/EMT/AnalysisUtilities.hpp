#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

#include <magic_enum/magic_enum.hpp>
#include <nlohmann/json.hpp>

#include <GridKit/CommonMath.hpp>
#include <GridKit/Model/EMT/SystemModelData.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace fs = ::std::filesystem;

    using Log = GridKit::Utilities::Logger;

    struct SwitchEvent
    {
      std::string element_id;
      bool        open;
    };

    struct SignalStep
    {
      std::string signal_id;
      double      value;
    };

    struct SystemEvent
    {
      double                                time;
      std::variant<SwitchEvent, SignalStep> action;
    };

    /**
     * @brief Data defined in JSON file for parameterized study
     */
    struct StudyData
    {
      std::map<std::string, std::map<std::string, double>> state;
      /// path to system model JSON file
      fs::path                                             system_model_file;
      /// path to model state JSON file, empty for none
      fs::path                                             state_file;
      /// monitor output time step size, or 0 for no intermediate monitoring
      double                                               dt_monitor;
      /// max time
      double                                               tmax;
      /// relative tolerance for the solver
      double                                               rel_tol;
      /// absolute tolerance for the solver
      double                                               abs_tol;
      /// Process-wide CommonMath smoothing scale
      double                                               mu{Math::DEFAULT_MU<double>};
      /// Study overrides for declared constant signals, by component path
      std::map<std::string, double>                        signal_values;
      /// fixed solver time step size, or 0 for adaptive stepping
      double                                               dt_fixed;
      /// maximum number of solver time steps, or 0 for the IDA default
      std::size_t                                          max_steps;
      /// IDA consistent initial condition calculation type
      AnalysisManager::Sundials::IdaConsistentICType       consistent_ic_type;
      /// set of system events
      std::vector<SystemEvent>                             events;
      /// path to output file
      fs::path                                             output_file;
      /// Optional CSV containing every DAE variable and derivative
      fs::path                                             state_output_file;
      /// path to reference file for validation
      fs::path                                             reference_file;
      /// Error tolerance (between output file and reference file)
      std::vector<double>                                  error_tol;
      /// Type of total error (relative or absolute)
      Testing::ErrorType                                   error_type;
      /// Smallest value at which to scale for relative error
      double                                               abs_err_threshold;
      /// Instance of model data
      SystemModelData<>                                    model_data;
    };

    using json = ::nlohmann::json;

    inline void from_json(const json& j, SystemEvent& event)
    {
      if (!j.is_object() || j.size() != 4 || !j.at("time").is_number())
        throw std::invalid_argument("An event requires time, type, and its two action fields");
      j.at("time").get_to(event.time);
      const auto type = j.at("type").get<std::string>();
      if (type == "switch")
      {
        if (!j.at("open").is_boolean())
          throw std::invalid_argument("A switch event requires a Boolean open value");
        event.action = SwitchEvent{j.at("element_id").get<std::string>(), j.at("open").get<bool>()};
      }
      else if (type == "signal_step")
      {
        if (!j.at("value").is_number() || !std::isfinite(j.at("value").get<double>()))
          throw std::invalid_argument("A signal_step value must be finite");
        event.action = SignalStep{j.at("signal_id").get<std::string>(), j.at("value").get<double>()};
      }
      else
      {
        throw std::invalid_argument("Unknown EMT event type: " + type);
      }
    }

    inline void validateEventTimes(const std::vector<SystemEvent>& events, double tmax)
    {
      if (!std::isfinite(tmax) || tmax < 0.0)
        throw std::invalid_argument("tmax must be finite and nonnegative");
      double previous = 0.0;
      for (const auto& event : events)
      {
        if (!std::isfinite(event.time) || event.time < previous || event.time > tmax)
          throw std::invalid_argument("Event times must be finite, ordered, and within [0, tmax]");
        previous = event.time;
      }
    }

    /** Configure before constructing models, including PWM's cached horizon. */
    template <typename RealT>
    inline void configureCommonMath(const StudyData& study)
    {
      Math::MU<RealT> = static_cast<RealT>(study.mu);
    }

    inline constexpr double DEFAULT_SOLVER_REL_TOL   = 1.0e-7;
    inline constexpr double DEFAULT_SOLVER_ABS_TOL   = 1.0e-9;
    inline constexpr double DEFAULT_VERIFICATION_TOL = 1.0e-4;

    /**
     * @brief JSON parser implemntation for `StudyData`
     */
    void from_json(const json& j, StudyData& c)
    {
      using namespace magic_enum;

      j.at("system_model_file").get_to(c.system_model_file);
      if (j.contains("state_file"))
      {
        j.at("state_file").get_to(c.state_file);
      }
      c.dt_monitor = j.value("dt_monitor", 0.0);
      j.at("tmax").get_to(c.tmax);
      c.rel_tol = j.value("rel_tol", DEFAULT_SOLVER_REL_TOL);
      c.abs_tol = j.value("abs_tol", DEFAULT_SOLVER_ABS_TOL);
      c.mu      = j.value("mu", Math::DEFAULT_MU<double>);
      if (!std::isfinite(c.mu) || c.mu <= 0.0)
      {
        throw std::invalid_argument("\"mu\" must be a positive finite number");
      }
      c.signal_values.clear();
      if (j.contains("signal_values"))
      {
        if (!j.at("signal_values").is_object())
        {
          throw std::invalid_argument("signal_values must be an object");
        }
        for (const auto& [name, value] : j.at("signal_values").items())
        {
          if (!value.is_number() || !std::isfinite(value.get<double>()))
          {
            throw std::invalid_argument("A signal_values entry must be a finite number");
          }
          c.signal_values.emplace(name, value.get<double>());
        }
      }
      c.dt_fixed           = j.value("dt_fixed", 0.0);
      c.max_steps          = j.value("max_steps", std::size_t{0});
      c.consistent_ic_type = AnalysisManager::Sundials::IdaConsistentICType::YA_YDP;
      if (j.contains("consistent_ic_type"))
      {
        const auto consistent_ic_type_str = j.at("consistent_ic_type").get<std::string>();
        if (consistent_ic_type_str == "y")
        {
          c.consistent_ic_type = AnalysisManager::Sundials::IdaConsistentICType::Y;
        }
        else if (consistent_ic_type_str == "ya_ydp")
        {
          c.consistent_ic_type = AnalysisManager::Sundials::IdaConsistentICType::YA_YDP;
        }
        else
        {
          Log::error() << "Invalid IDA consistent initial condition type \""
                       << consistent_ic_type_str << "\"; "
                       << "must be either \"y\" or \"ya_ydp\"";
        }
      }

      c.events.clear();
      if (j.contains("events"))
        j.at("events").get_to(c.events);
      validateEventTimes(c.events, c.tmax);

      if (j.contains("output_file"))
      {
        j.at("output_file").get_to(c.output_file);
      }
      if (j.contains("state_output_file"))
      {
        j.at("state_output_file").get_to(c.state_output_file);
      }

      if (j.contains("reference_file"))
      {
        j.at("reference_file").get_to(c.reference_file);
      }

      if (j.contains("error_tolerance"))
      {
        auto& tolj = j.at("error_tolerance");
        if (tolj.is_array())
        {
          tolj.get_to(c.error_tol);
        }
        else
        {
          tolj.get_to(c.error_tol.emplace_back());
        }
      }
      else
      {
        c.error_tol.push_back(DEFAULT_VERIFICATION_TOL);
      }

      using ErrorType = Testing::ErrorType;
      if (j.contains("error_type"))
      {
        auto type_str  = j.at("error_type").get<std::string>();
        auto type_wrap = enum_cast<ErrorType>(type_str, case_insensitive);
        if (!type_wrap.has_value())
        {
          Log::error() << "Invalid error type \"" << type_str << "\"; "
                       << "must be either \"relative\" or \"absolute\"";
        }
        c.error_type = type_wrap.value();
      }
      else
      {
        c.error_type = ErrorType::RELATIVE;
      }

      c.abs_err_threshold = j.value("abs_err_threshold", Testing::DEFAULT_ABS_ERROR_THRESHOLD);
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
      if (!data.state_file.empty())
      {
        if (!data.state_file.is_absolute())
        {
          data.state_file = loc / data.state_file;
        }
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

      // Override only explicitly declared constants, never component outputs.
      for (const auto& [path, value] : data.signal_values)
      {
        ContainerData<double, size_t>* scope     = &data.model_data;
        auto                           remaining = path;
        while (remaining.find('.') != std::string::npos)
        {
          const auto dot   = remaining.find('.');
          const auto name  = remaining.substr(0, dot);
          auto       child = std::find_if(scope->container.begin(), scope->container.end(), [&](const auto& candidate)
                                    { return candidate.id == name; });
          if (child == scope->container.end())
          {
            throw std::invalid_argument("Unknown constant signal: " + path);
          }
          scope = &*child;
          remaining.erase(0, dot + 1);
        }
        auto signal = std::find_if(scope->signal.begin(), scope->signal.end(), [&](const auto& candidate)
                                   { return candidate.id == remaining; });
        if (signal == scope->signal.end() || !signal->value.has_value())
        {
          throw std::invalid_argument("signal_values requires a declared constant: " + path);
        }
        signal->value = value;
      }

      if (!data.state_file.empty())
      {
        std::ifstream state_stream(data.state_file);
        if (!state_stream)
          throw std::invalid_argument("Cannot open state file: " + data.state_file.string());
        const auto state = json::parse(state_stream);
        if (!state.is_object())
          throw std::invalid_argument("State must be a JSON object");
        for (const auto* section : {"buses", "devices"})
        {
          if (!state.contains(section) || state.at(section).is_null())
            continue;
          if (!state.at(section).is_object())
            throw std::invalid_argument(std::string("State ") + section + " must be an object");
          for (const auto& [path, outputs] : state.at(section).items())
          {
            if (outputs.is_null())
              continue;
            if (!outputs.is_object())
              throw std::invalid_argument("Component state must be an object: " + path);
            for (const auto& [name, value] : outputs.items())
            {
              if (value.is_null() || name == "injections")
                continue;
              if (name == "open" && value.is_boolean())
                data.state[path][name] = value.get<bool>() ? 1.0 : 0.0;
              else if (value.is_number() && std::isfinite(value.get<double>()))
                data.state[path][name] = value.get<double>();
              else
                throw std::invalid_argument("Initial output must be finite: " + path + "." + name);
            }
          }
        }
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

      const auto& out_file = study_data.output_file;
      const auto& ref_file = study_data.reference_file;
      if (!out_file.empty() && !ref_file.empty())
      {
        auto errorSet = Testing::compareCSV(out_file,
                                            ref_file,
                                            study_data.error_type,
                                            study_data.abs_err_threshold);

        // Print the errors
        if (print_results)
        {
          errorSet->display();
        }

        // Check against specified tolerance
        if (study_data.error_tol.size() > 1)
        {
          status *= study_data.error_tol.size() == errorSet->var_errors.size();
          for (std::size_t i = 0; i < study_data.error_tol.size(); ++i)
          {
            status *= errorSet->var_errors[i].max_value < study_data.error_tol[i];
          }
        }
        else
        {
          status *= errorSet->total_error.max_value < study_data.error_tol[0];
        }

        if (print_results)
        {
          status.report();
        }
      }
      return status;
    }

  } // namespace EMT
} // namespace GridKit
