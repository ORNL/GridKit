#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <filesystem>
#include <format>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include <magic_enum/magic_enum.hpp>
#include <nlohmann/json.hpp>

#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>
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
     * @brief Data defined in JSON file for parameterized study
     */
    struct StudyData
    {
      /// path to system model JSON file
      fs::path                                       system_model_file;
      /// monitor output time step size, or 0 for no intermediate monitoring
      double                                         dt_monitor;
      /// max time
      double                                         tmax;
      /// IDA solver options
      AnalysisManager::Sundials::IdaOptions<double>  ida;
      /// IDA consistent initial condition calculation type
      AnalysisManager::Sundials::IdaConsistentICType consistent_ic_type;
      /// bus where the study's first bus fault is applied
      std::optional<std::size_t>                     fault_bus;
      /// set of system events
      std::vector<SystemEvent>                       events;
      /// path to output file
      fs::path                                       output_file;
      /// path to per-step integrator trace file, or empty for no trace
      fs::path                                       step_trace_file;
      /// path to reference file for validation
      fs::path                                       reference_file;
      /// Error tolerance (between output file and reference file)
      std::vector<double>                            error_tol;
      /// Type of total error (relative or absolute)
      Testing::ErrorType                             error_type;
      /// Smallest value at which to scale for relative error
      double                                         abs_err_threshold;
      /// Instance of model data
      SystemModelData<>                              model_data;
    };

    using json = ::nlohmann::json;
    using Log  = ::GridKit::Utilities::Logger;

    inline constexpr double DEFAULT_SOLVER_REL_TOL   = 1.0e-7;
    inline constexpr double DEFAULT_SOLVER_ABS_TOL   = 1.0e-9;
    inline constexpr double DEFAULT_VERIFICATION_TOL = 1.0e-4;

    inline constexpr std::array<std::string_view, 23> IDA_OPTION_KEYS = {
        "rel_tol",
        "abs_tol",
        "fixed_step",
        "init_step",
        "min_step",
        "max_step",
        "max_order",
        "max_num_steps",
        "max_err_test_fails",
        "suppress_alg",
        "max_nonlin_iters",
        "max_conv_fails",
        "nonlin_conv_coef",
        "max_num_steps_ic",
        "max_num_jacs_ic",
        "max_num_iters_ic",
        "max_backs_ic",
        "line_search_off_ic",
        "nonlin_conv_coef_ic",
        "step_tolerance_ic",
        "linear_solution_scaling",
        "delta_cj_lsetup",
        "klu_ordering"};

    template <class T>
    void getOptional(const json& j, std::string_view key, std::optional<T>& value)
    {
      if (j.contains(key))
      {
        value = j.at(key).get<T>();
      }
    }

    void parseKluOrdering(const json&                                            value,
                          std::optional<AnalysisManager::Sundials::KluOrdering>& ordering)
    {
      using KluOrdering = AnalysisManager::Sundials::KluOrdering;
      const auto name   = value.get<std::string>();
      const auto parsed = magic_enum::enum_cast<KluOrdering>(
          name,
          magic_enum::case_insensitive);
      if (!parsed.has_value())
      {
        throw std::invalid_argument("klu_ordering must be amd, colamd, or natural");
      }
      ordering = *parsed;
    }

    void parseIdaOptions(const json&                                    j,
                         AnalysisManager::Sundials::IdaOptions<double>& options)
    {
      if (!j.is_object())
      {
        throw std::invalid_argument("ida must be a JSON object");
      }

      for (auto entry = j.begin(); entry != j.end(); ++entry)
      {
        if (std::find(IDA_OPTION_KEYS.begin(), IDA_OPTION_KEYS.end(), entry.key())
            == IDA_OPTION_KEYS.end())
        {
          throw std::invalid_argument("Unknown IDA option: " + entry.key());
        }
      }

      options.rel_tol = j.value("rel_tol", options.rel_tol);
      options.abs_tol = j.value("abs_tol", options.abs_tol);

      getOptional(j, "fixed_step", options.fixed_step);
      getOptional(j, "init_step", options.init_step);
      getOptional(j, "min_step", options.min_step);
      getOptional(j, "max_step", options.max_step);
      getOptional(j, "max_order", options.max_order);
      getOptional(j, "max_num_steps", options.max_num_steps);
      getOptional(j, "max_err_test_fails", options.max_err_test_fails);
      getOptional(j, "suppress_alg", options.suppress_alg);
      getOptional(j, "max_nonlin_iters", options.max_nonlin_iters);
      getOptional(j, "max_conv_fails", options.max_conv_fails);
      getOptional(j, "nonlin_conv_coef", options.nonlin_conv_coef);
      getOptional(j, "max_num_steps_ic", options.max_num_steps_ic);
      getOptional(j, "max_num_jacs_ic", options.max_num_jacs_ic);
      getOptional(j, "max_num_iters_ic", options.max_num_iters_ic);
      getOptional(j, "max_backs_ic", options.max_backs_ic);
      getOptional(j, "line_search_off_ic", options.line_search_off_ic);
      getOptional(j, "nonlin_conv_coef_ic", options.nonlin_conv_coef_ic);
      getOptional(j, "step_tolerance_ic", options.step_tolerance_ic);
      getOptional(j, "linear_solution_scaling", options.linear_solution_scaling);
      getOptional(j, "delta_cj_lsetup", options.delta_cj_lsetup);

      if (j.contains("klu_ordering"))
      {
        parseKluOrdering(j.at("klu_ordering"), options.klu_ordering);
      }
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

      c.ida.rel_tol      = DEFAULT_SOLVER_REL_TOL;
      c.ida.abs_tol      = DEFAULT_SOLVER_ABS_TOL;
      c.ida.klu_ordering = AnalysisManager::Sundials::KluOrdering::AMD;

      const bool has_legacy_ida = j.contains("rel_tol")
                                  || j.contains("abs_tol")
                                  || j.contains("dt_fixed")
                                  || j.contains("max_steps")
                                  || j.contains("klu_ordering");
      if (has_legacy_ida && j.contains("ida"))
      {
        throw std::invalid_argument(
            "Specify IDA options either at the top level or in ida, not both");
      }

      if (j.contains("ida"))
      {
        parseIdaOptions(j.at("ida"), c.ida);
      }
      else
      {
        c.ida.rel_tol = j.value("rel_tol", c.ida.rel_tol);
        c.ida.abs_tol = j.value("abs_tol", c.ida.abs_tol);

        if (j.contains("dt_fixed"))
        {
          const auto fixed_step = j.at("dt_fixed").get<double>();
          if (fixed_step != 0.0)
          {
            c.ida.fixed_step = fixed_step;
          }
        }
        if (j.contains("max_steps"))
        {
          const auto max_steps = j.at("max_steps").get<std::size_t>();
          if (max_steps != 0)
          {
            c.ida.max_num_steps = static_cast<long int>(max_steps);
          }
        }
        if (j.contains("klu_ordering"))
        {
          parseKluOrdering(j.at("klu_ordering"), c.ida.klu_ordering);
        }
      }

      c.consistent_ic_type = AnalysisManager::Sundials::IdaConsistentICType::YA_YDP;
      if (j.contains("fault_bus"))
      {
        c.fault_bus = j.at("fault_bus").get<std::size_t>();
      }
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

      if (j.contains("step_trace_file"))
      {
        j.at("step_trace_file").get_to(c.step_trace_file);
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
      if (!data.reference_file.empty())
      {
        if (!data.reference_file.is_absolute())
        {
          data.reference_file = loc / data.reference_file;
        }
      }
      if (!data.step_trace_file.empty() && !data.step_trace_file.is_absolute())
      {
        data.step_trace_file = loc / data.step_trace_file;
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

    /**
     * @brief Write the integrator's per-step trace as CSV
     */
    void writeStepTrace(
        const fs::path&                                              file_path,
        const std::vector<AnalysisManager::Sundials::IdaStepRecord>& trace)
    {
      auto out = std::ofstream(file_path);
      if (!out)
      {
        throw std::runtime_error("Failed to open step trace file: " + file_path.string());
      }

      out << "segment,t,h,h_next,order,order_next,nsteps,nres,njac,netf,nni,nncf\n";
      out << std::setprecision(12);
      for (const auto& record : trace)
      {
        out << record.segment << ','
            << record.t << ','
            << record.h << ','
            << record.h_next << ','
            << record.order << ','
            << record.order_next << ','
            << record.num_steps << ','
            << record.num_residual_evals << ','
            << record.num_jacobian_evals << ','
            << record.num_error_test_fails << ','
            << record.num_nonlinear_iters << ','
            << record.num_nonlinear_convergence_fails << '\n';
      }

      out.flush();
      if (!out)
      {
        throw std::runtime_error("Failed to write step trace file: " + file_path.string());
      }
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

  } // namespace PhasorDynamics
} // namespace GridKit
