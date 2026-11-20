/**
 * @file AnalysisUtilities.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Study files of the EMT transient simulation.
 */

#pragma once

#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <format>
#include <fstream>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include <magic_enum/magic_enum.hpp>
#include <nlohmann/json.hpp>

#include <GridKit/Model/EMT/SystemModelData.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>
#include <GridKit/Testing/AggregateErrors.hpp>
#include <GridKit/Testing/CSV.hpp>
#include <GridKit/Testing/Tokenizer.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace EMT
  {
    namespace fs = ::std::filesystem;

    using json = ::nlohmann::json;
    using Log  = ::GridKit::Utilities::Logger;

    struct StudyData
    {
      fs::path                                       system_model_file;
      fs::path                                       state_file; ///< Empty for a zero state
      double                                         dt_monitor; ///< Zero for no intermediate output
      double                                         tmax;
      double                                         rel_tol;
      double                                         abs_tol;   ///< Zero for the tolerances of the models
      double                                         dt_fixed;  ///< Zero for adaptive steps
      std::size_t                                    max_steps; ///< Steps between outputs
      int                                            max_order;
      bool                                           suppress_algebraic_errors; ///< Leave the algebraic variables out of the error test
      AnalysisManager::Sundials::IdaConsistentICType consistent_ic_type;
      std::vector<json>                              events; ///< State values to set, each at its `time`
      fs::path                                       output_file;
      fs::path                                       reference_file; ///< Empty for no comparison
      double                                         error_tol;      ///< Bound on the total relative error against the reference
      ::GridKit::Model::VariableMonitorLayout        output_layout;
      SystemModelData<>                              model_data;
      json                                           state; ///< Instantaneous values at t0
    };

    inline constexpr double DEFAULT_SOLVER_REL_TOL   = 1.0e-2;
    inline constexpr double DEFAULT_VERIFICATION_TOL = 1.0e-4;

    /// Largest difference between a reference time and the time of the monitor row it is compared with
    inline constexpr double REFERENCE_TIME_TOL = 1.0e-9;

    /// IDA counts steps in `long`, so no limit
    inline constexpr std::size_t UNLIMITED_STEPS = std::numeric_limits<long>::max();

    inline void from_json(const json& j, StudyData& c)
    {
      using namespace magic_enum;

      j.at("system_model_file").get_to(c.system_model_file);
      c.state_file = j.value("state_file", fs::path{});
      c.dt_monitor = j.value("dt_monitor", 0.0);
      j.at("tmax").get_to(c.tmax);
      c.rel_tol                   = j.value("rel_tol", DEFAULT_SOLVER_REL_TOL);
      c.abs_tol                   = j.value("abs_tol", 0.0);
      c.dt_fixed                  = j.value("dt_fixed", 0.0);
      c.max_steps                 = j.value("max_steps", UNLIMITED_STEPS);
      c.max_order                 = j.value("max_order", 2);
      c.suppress_algebraic_errors = j.value("suppress_algebraic_errors", true);
      c.consistent_ic_type        = AnalysisManager::Sundials::IdaConsistentICType::YA_YDP;

      c.events = j.value("events", std::vector<json>{});

      c.output_file    = j.value("output_file", fs::path{});
      c.reference_file = j.value("reference_file", fs::path{});
      c.error_tol      = j.value("error_tolerance", DEFAULT_VERIFICATION_TOL);

      using OutputLayout = ::GridKit::Model::VariableMonitorLayout;
      c.output_layout    = OutputLayout::WIDE;
      if (j.contains("output_layout"))
      {
        const auto layout_str  = j.at("output_layout").get<std::string>();
        const auto layout_wrap = enum_cast<OutputLayout>(layout_str, case_insensitive);
        if (!layout_wrap.has_value())
        {
          throw std::invalid_argument(std::format("Invalid output layout \"{}\"", layout_str));
        }
        c.output_layout = layout_wrap.value();
      }
    }

    /// The case, state and reference paths are relative to the study file
    inline StudyData parseStudyData(const fs::path& file_path)
    {
      std::ifstream stream(file_path);
      if (!stream)
      {
        throw std::runtime_error("Could not open file " + file_path.string());
      }
      auto data = json::parse(stream).get<StudyData>();

      const auto loc = file_path.parent_path();
      if (!data.system_model_file.is_absolute())
      {
        data.system_model_file = loc / data.system_model_file;
      }
      data.model_data = parseSystemModelData(data.system_model_file);

      data.state = json::object();
      if (!data.state_file.empty())
      {
        if (!data.state_file.is_absolute())
        {
          data.state_file = loc / data.state_file;
        }
        std::ifstream state(data.state_file);
        if (!state)
        {
          throw std::runtime_error("Could not open file " + data.state_file.string());
        }
        data.state = json::parse(state);
      }

      if (!data.reference_file.empty() && !data.reference_file.is_absolute())
      {
        data.reference_file = loc / data.reference_file;
      }

      if (!data.output_file.empty())
      {
        data.model_data.monitor_sink.emplace_back(::GridKit::Model::VariableMonitorFormat::CSV,
                                                  data.output_file.string(),
                                                  ",",
                                                  data.output_layout);
      }
      return data;
    }

    inline void checkCommandLine(int argc, const std::string& appName)
    {
      if (argc < 2)
      {
        Log::error() << "No input file provided\n"
                     << std::format(
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

    /**
     * @brief Compare the monitor output with the reference file.
     *
     * The reference holds monitored columns, found by name, at monitor times.
     * It is compared with the first monitor row at each of its times, since an
     * event or a current zero adds rows at its own time.
     *
     * @param[in] study Study with a wide CSV output file.
     * @return Whether the total relative error is below the tolerance of the
     *         study, or a pass when the study has no reference.
     *
     * @throws std::runtime_error if a reference column or time is not in the output.
     */
    inline Testing::TestStatus checkErrors(const StudyData& study)
    {
      auto status = Testing::TestStatus{"monitor file vs reference file"};
      if (study.output_file.empty() || study.reference_file.empty())
      {
        return status;
      }

      auto output    = Testing::checkOpenFile(study.output_file.string());
      auto reference = Testing::checkOpenFile(study.reference_file.string());

      std::string line;
      std::getline(output, line);
      const auto monitored = Testing::Tokenizer<>(line, ',')();
      std::getline(reference, line);
      const auto labels = Testing::Tokenizer<>(line, ',')();

      // Column of the monitor output for each reference column, the time first
      std::vector<std::size_t> columns{0};
      for (const auto& label : std::span{labels}.subspan(1))
      {
        const auto found = std::ranges::find(monitored, label);
        if (found == monitored.end())
        {
          throw std::runtime_error(std::format("Reference column {} is not monitored", label));
        }
        columns.push_back(static_cast<std::size_t>(found - monitored.begin()));
      }

      auto errors = Testing::makeErrorSet(Testing::ErrorType::RELATIVE, std::span{labels}.subspan(1));
      while (std::getline(reference, line))
      {
        const Testing::OutputAtTime expected(Testing::Tokenizer<double>(line, ',')());

        std::vector<double> row;
        do
        {
          if (!std::getline(output, line))
          {
            throw std::runtime_error(std::format("No monitor output at reference time {}", expected.t));
          }
          row = Testing::Tokenizer<double>(line, ',')();
        } while (row.front() < expected.t - REFERENCE_TIME_TOL);
        if (row.front() > expected.t + REFERENCE_TIME_TOL)
        {
          throw std::runtime_error(std::format("No monitor output at reference time {}", expected.t));
        }

        std::vector<double> values;
        for (const auto column : columns)
        {
          values.push_back(row[column]);
        }
        errors->push(Testing::OutputAtTime(values), expected);
      }
      errors->wrap();
      errors->display();

      status *= errors->total_error.max_value < study.error_tol;
      status.report();
      return status;
    }
  } // namespace EMT
} // namespace GridKit
