#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <variant>
#include <vector>

#include <magic_enum/magic_enum.hpp>
#include <nlohmann/json.hpp>

#include <GridKit/CommonMath.hpp>
#include <GridKit/Model/EMT/JsonValidation.hpp>
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
     * @brief Reflected-current data defining a constant or harmonic line prehistory
     */
    struct LineHistory
    {
      double                           omega;
      std::array<ABCVector<double>, 2> value, derivative;
    };

    /**
     * @brief Data defined in JSON file for parameterized study
     */
    struct StudyData
    {
      std::map<std::string, LineHistory>                   history;
      std::map<std::string, std::map<std::string, double>> state;
      /// Angular frequency of a balanced initial operating point; zero uses explicit output values.
      double                                               initial_omega{0.0};
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
      bool                                                 scaled_abs_tol{false};
      /// Process-wide CommonMath smoothing scale
      double                                               mu{Math::DEFAULT_MU<double>};
      /// Study overrides for declared constant signals, by component path
      std::map<std::string, double>                        signal_values;
      /// fixed solver time step size, or 0 for adaptive stepping
      double                                               dt_fixed;
      /// maximum number of solver time steps, or 0 for the IDA default
      std::size_t                                          max_steps;
      /// Maximum BDF order for adaptive integration
      int                                                  max_order{5};
      /// IDA consistent initial condition calculation type
      AnalysisManager::Sundials::IdaConsistentICType       consistent_ic_type;
      /// set of system events
      std::vector<SystemEvent>                             events;
      /// path to output file
      fs::path                                             output_file;
      /// Optional CSV containing every DAE variable and derivative
      fs::path                                             state_output_file;
      /// Optional accepted internal step log, independent of monitor cadence
      fs::path                                             step_output_file;
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

      validateJsonFields(j, "EMT study", {"system_model_file", "state_file", "dt_monitor", "tmax", "rel_tol", "abs_tol", "scaled_abs_tol", "mu", "signal_values", "dt_fixed", "max_steps", "max_order", "consistent_ic_type", "events", "output_file", "state_output_file", "step_output_file", "reference_file", "error_tolerance", "error_type", "abs_err_threshold"});
      const auto real = [&j](const char* key, double fallback)
      { return j.contains(key) ? parseFiniteReal<double>(j.at(key), key) : fallback; };

      j.at("system_model_file").get_to(c.system_model_file);
      if (j.contains("state_file"))
      {
        j.at("state_file").get_to(c.state_file);
      }
      c.dt_monitor     = real("dt_monitor", 0.0);
      c.tmax           = parseFiniteReal<double>(j.at("tmax"), "tmax");
      c.rel_tol        = real("rel_tol", DEFAULT_SOLVER_REL_TOL);
      c.abs_tol        = real("abs_tol", DEFAULT_SOLVER_ABS_TOL);
      c.scaled_abs_tol = false;
      if (j.contains("scaled_abs_tol"))
      {
        if (!j.at("scaled_abs_tol").is_boolean())
          throw std::invalid_argument("scaled_abs_tol must be a boolean");
        c.scaled_abs_tol = j.at("scaled_abs_tol").get<bool>();
      }
      c.mu = real("mu", Math::DEFAULT_MU<double>);
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
          c.signal_values.emplace(name, parseFiniteReal<double>(value, "signal_values entry \"" + name + "\""));
        }
      }
      c.dt_fixed = real("dt_fixed", 0.0);
      if (c.dt_monitor < 0 || c.dt_fixed < 0 || c.rel_tol <= 0 || c.abs_tol <= 0)
        throw std::invalid_argument("EMT study requires nonnegative time steps and positive solver tolerances");
      c.max_steps = 0;
      if (j.contains("max_steps"))
      {
        const auto& steps = j.at("max_steps");
        if (!steps.is_number_integer()
            || (!steps.is_number_unsigned() && steps.get<int64_t>() < 0)
            || steps.get<uint64_t>() > static_cast<uint64_t>(std::numeric_limits<long int>::max()))
          throw std::invalid_argument("max_steps requires a nonnegative integer within the solver's index range");
        c.max_steps = steps.get<size_t>();
      }
      c.max_order = 5;
      if (j.contains("max_order"))
      {
        const auto& order = j.at("max_order");
        if (!order.is_number_integer() || order < 1 || order > 5)
          throw std::invalid_argument("max_order requires an integer between 1 and 5");
        c.max_order = order.get<int>();
      }
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
          throw std::invalid_argument("Invalid consistent_ic_type \"" + consistent_ic_type_str
                                      + "\"; expected \"y\" or \"ya_ydp\"");
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

      if (j.contains("step_output_file"))
        j.at("step_output_file").get_to(c.step_output_file);

      if (j.contains("reference_file"))
      {
        j.at("reference_file").get_to(c.reference_file);
      }

      c.error_tol.clear();
      if (j.contains("error_tolerance"))
      {
        auto& tolj = j.at("error_tolerance");
        if (tolj.is_array())
        {
          if (tolj.empty())
            throw std::invalid_argument("error_tolerance must not be empty");
          for (const auto& value : tolj)
            c.error_tol.push_back(parseFiniteReal<double>(value, "error_tolerance"));
        }
        else
        {
          c.error_tol.push_back(parseFiniteReal<double>(tolj, "error_tolerance"));
        }
      }
      else
      {
        c.error_tol.push_back(DEFAULT_VERIFICATION_TOL);
      }
      if (std::any_of(c.error_tol.begin(), c.error_tol.end(), [](double tolerance)
                      { return tolerance <= 0; }))
        throw std::invalid_argument("error_tolerance must be positive");

      using ErrorType = Testing::ErrorType;
      if (j.contains("error_type"))
      {
        auto type_str  = j.at("error_type").get<std::string>();
        auto type_wrap = enum_cast<ErrorType>(type_str, case_insensitive);
        if (!type_wrap.has_value())
        {
          throw std::invalid_argument("Invalid error_type \"" + type_str + "\"; expected \"relative\" or \"absolute\"");
        }
        c.error_type = type_wrap.value();
      }
      else
      {
        c.error_type = ErrorType::RELATIVE;
      }

      c.abs_err_threshold = real("abs_err_threshold", Testing::DEFAULT_ABS_ERROR_THRESHOLD);
      if (c.abs_err_threshold < 0)
        throw std::invalid_argument("abs_err_threshold must be nonnegative");
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

    /// Validate all named records before discarding null values that request model defaults.
    inline std::map<std::string, std::map<std::string, double>>
    parseInitialState(const json& state, const ContainerData<double, size_t>& model)
    {
      validateJsonFields(state, "State", {"header", "buses", "devices", "history"});
      if (state.contains("header") && !state.at("header").is_null())
      {
        const auto& header = state.at("header");
        validateJsonFields(header, "State header", {"version", "time", "omega", "created", "description"});
        if (header.contains("omega") && !header.at("omega").is_null() && parseFiniteReal<double>(header.at("omega"), "Initial angular frequency") <= 0.0)
          throw std::invalid_argument("Initial angular frequency must be positive");
        if (header.contains("version") && !header.at("version").is_null())
        {
          const auto& version = header.at("version");
          if (!version.is_number_integer()
              || (!version.is_number_unsigned() && version.get<int64_t>() < 0)
              || version.get<uint64_t>() > std::numeric_limits<unsigned int>::max())
            throw std::invalid_argument("State header version requires a nonnegative integer within the version range");
        }
        if (header.contains("time") && !header.at("time").is_null()
            && parseFiniteReal<double>(header.at("time"), "State header time") != 0)
          throw std::invalid_argument("The EMT application requires an initial state at time zero");
        for (const auto* name : {"created", "description"})
          if (header.contains(name) && !header.at(name).is_null() && !header.at(name).is_string())
            throw std::invalid_argument(std::string("State header ") + name + " must be a string");
      }

      std::map<std::string, std::set<std::string>> allowed;
      std::set<std::string>                        buses;
      const auto                                   collect = [&](auto&& self, const auto& scope, const std::string& prefix) -> void
      {
        const auto add = [&](const auto& devices)
        {
          using Data    = typename std::decay_t<decltype(devices)>::value_type;
          using Outputs = typename Data::Outputs;
          for (const auto& device : devices)
          {
            const auto path  = prefix + device.id;
            auto&      names = allowed[path];
            for (const auto output : magic_enum::enum_values<Outputs>())
              if (output != Outputs::SIZE)
                names.emplace(magic_enum::enum_name(output));
            if constexpr (std::is_same_v<Outputs, SwitchOutputs>)
              names.insert("open");
            if constexpr (std::is_same_v<Outputs, TransformerOutputs>)
              names = {"i12a", "i12b", "i12c", "psi1a", "psi1b", "psi1c", "psi2a", "psi2b", "psi2c"};
            if constexpr (std::is_same_v<Outputs, BusOutputs>)
              buses.insert(path);
          }
        };
        std::apply([&](const auto&... devices)
                   { (add(devices), ...); },
                   std::tie(scope.bus, scope.loadz, scope.voltage_source, scope.dependent_voltage_source, scope.filter, scope.machine, scope.regfma, scope.line_lumped, scope.line_distributed, scope.sw, scope.transformer, scope.inner_current_control, scope.park, scope.pll, scope.pwm, scope.converter, scope.ieeest, scope.gastpti, scope.gov, scope.sexs_pti, scope.exciter));
        for (const auto& child : scope.container)
          self(self, child, prefix + child.id + ".");
      };
      collect(collect, model, "");

      std::map<std::string, std::map<std::string, double>> result;
      std::set<std::string>                                seen;
      for (const auto* section : {"buses", "devices"})
      {
        if (!state.contains(section) || state.at(section).is_null())
          continue;
        if (!state.at(section).is_object())
          throw std::invalid_argument(std::string("State ") + section + " must be an object");
        for (const auto& [path, outputs] : state.at(section).items())
        {
          const auto record = allowed.find(path);
          if (record == allowed.end())
            throw std::invalid_argument("Unknown initial state component: " + path);
          if (std::string_view(section) == "buses" && !buses.contains(path))
            throw std::invalid_argument("Initial bus state requires a Bus component: " + path);
          if (!seen.insert(path).second)
            throw std::invalid_argument("Duplicate initial state component: " + path);
          if (outputs.is_null())
            continue;
          if (!outputs.is_object())
            throw std::invalid_argument("Component state must be an object: " + path);
          for (const auto& [name, value] : outputs.items())
          {
            if (!record->second.contains(name))
              throw std::invalid_argument("Unknown initial output: " + path + "." + name);
            if (value.is_null())
              continue;
            if (name == "open")
            {
              if (!value.is_boolean())
                throw std::invalid_argument("Initial switch open must be Boolean: " + path);
              result[path][name] = value.get<bool>() ? 1.0 : 0.0;
            }
            else
              result[path][name] = parseFiniteReal<double>(value, "Initial output " + path + "." + name);
          }
        }
      }
      return result;
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
        data.state       = parseInitialState(state, data.model_data);
        if (state.contains("header") && state.at("header").is_object()
            && state.at("header").contains("omega") && !state.at("header").at("omega").is_null())
          data.initial_omega = state.at("header").at("omega").get<double>();
        if (state.contains("history"))
        {
          if (!state.at("history").is_object())
            throw std::invalid_argument("State history must be an object");
          for (const auto& [path, entry] : state.at("history").items())
          {
            validateJsonFields(entry, "Line history " + path, {"omega", "i_ref1", "i_ref2", "d_i_ref1", "d_i_ref2"});
            LineHistory history;
            history.omega = parseFiniteReal<double>(entry.at("omega"), "History omega");
            for (size_t e = 0; e < 2; ++e)
            {
              const auto  name       = "i_ref" + std::to_string(e + 1);
              const auto& value      = entry.at(name);
              const auto& derivative = entry.at("d_" + name);
              if (!value.is_array() || value.size() != 3 || !derivative.is_array() || derivative.size() != 3)
                throw std::invalid_argument("Line history requires three values and derivatives per terminal");
              for (size_t p = 0; p < 3; ++p)
              {
                history.value[e][p]      = parseFiniteReal<double>(value[p], "History value");
                history.derivative[e][p] = parseFiniteReal<double>(derivative[p], "History derivative");
              }
            }
            data.history.emplace(path, history);
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
