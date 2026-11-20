/**
 * @file AnalysisUtilities.hpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Study files of the EMT transient simulation.
 */

#pragma once

#include <cstddef>
#include <filesystem>
#include <format>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <magic_enum/magic_enum.hpp>
#include <nlohmann/json.hpp>

#include <GridKit/Model/EMT/SystemModelData.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>
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
      ::GridKit::Model::VariableMonitorLayout        output_layout;
      SystemModelData<>                              model_data;
      json                                           state; ///< Instantaneous values at t0
    };

    inline constexpr double DEFAULT_SOLVER_REL_TOL = 1.0e-2;

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

      c.output_file = j.value("output_file", fs::path{});

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

    /// The case and state paths are relative to the study file
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
  } // namespace EMT
} // namespace GridKit
