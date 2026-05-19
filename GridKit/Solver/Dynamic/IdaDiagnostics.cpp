#include <fstream>
#include <stdexcept>
#include <string>
#include <utility>

#include <nlohmann/json.hpp>

#include <GridKit/Solver/Dynamic/IdaDiagnostics.hpp>

namespace AnalysisManager
{
  namespace Sundials
  {
    namespace
    {
      std::string logLevelName(IdaLogLevel level)
      {
        return level == IdaLogLevel::Error ? "error" : "warning";
      }

      nlohmann::json idaStatsJson(const IdaStats& stats)
      {
        return nlohmann::json{
            {"sundials",
             {{"version", stats.sundials_version_},
              {"logging_level", stats.sundials_logging_level_}}},
            {"integrator",
             {{"steps", stats.num_steps_},
              {"residual_evals", stats.num_residual_evals_},
              {"linear_solver_setups", stats.num_linear_solver_setups_},
              {"error_test_failures", stats.num_error_test_fails_},
              {"backtrack_operations", stats.num_backtrack_operations_}}},
            {"nonlinear_solver",
             {{"iterations", stats.num_nonlinear_iters_},
              {"convergence_failures", stats.num_nonlinear_convergence_fails_},
              {"step_solve_failures", stats.num_nonlinear_step_fails_}}},
            {"linear_solver",
             {{"jacobian_evals", stats.num_jacobian_evals_},
              {"last_jacobian_eval_step", stats.num_jacobian_eval_steps_},
              {"jacobian_time", stats.jacobian_time_},
              {"jacobian_cj", stats.jacobian_cj_},
              {"iterations", stats.num_linear_iters_},
              {"convergence_failures", stats.num_linear_convergence_fails_},
              {"residual_evals", stats.num_linear_residual_evals_},
              {"preconditioner_evals", stats.num_preconditioner_evals_},
              {"preconditioner_solves", stats.num_preconditioner_solves_},
              {"jtimes_setup_evals", stats.num_jtimes_setup_evals_},
              {"jtimes_evals", stats.num_jtimes_evals_},
              {"last_flag", stats.last_linear_flag_},
              {"last_flag_name", stats.last_linear_flag_name_}}},
            {"final_state",
             {{"last_order", stats.last_order_},
              {"current_order", stats.current_order_},
              {"actual_initial_step", stats.actual_initial_step_},
              {"last_step", stats.last_step_},
              {"current_step", stats.current_step_},
              {"current_time", stats.current_time_},
              {"current_cj", stats.current_cj_}}}};
      }

      nlohmann::json idaStatsSegmentJson(const IdaStatsSegment& segment)
      {
        auto json = idaStatsJson(segment.stats);
        json.erase("sundials");
        json["start_time"]   = segment.start_time;
        json["end_time"]     = segment.end_time;
        json["output_steps"] = segment.output_steps;
        return json;
      }
    } // namespace

    IdaStatsRecorder::IdaStatsRecorder(bool enabled)
      : enabled_(enabled)
    {
    }

    bool IdaStatsRecorder::enabled() const
    {
      return enabled_;
    }

    void IdaStatsRecorder::recordSegment(const IdaStats& start_stats,
                                         const IdaStats& end_stats,
                                         double          start_time,
                                         double          end_time,
                                         int             output_steps)
    {
      auto stats  = idaStatsDelta(end_stats, start_stats);
      summary_   += stats;
      segments_.push_back({start_time, end_time, output_steps, std::move(stats)});
    }

    IdaStatsReport IdaStatsRecorder::report(std::optional<IdaLogOptions> log) const
    {
      return {summary_, segments_, std::move(log)};
    }

    IdaStats idaStatsDelta(const IdaStats& end_stats, const IdaStats& start_stats)
    {
      IdaStats delta = end_stats;

      delta.num_steps_                       = end_stats.num_steps_ - start_stats.num_steps_;
      delta.num_residual_evals_              = end_stats.num_residual_evals_ - start_stats.num_residual_evals_;
      delta.num_linear_solver_setups_        = end_stats.num_linear_solver_setups_ - start_stats.num_linear_solver_setups_;
      delta.num_error_test_fails_            = end_stats.num_error_test_fails_ - start_stats.num_error_test_fails_;
      delta.num_backtrack_operations_        = end_stats.num_backtrack_operations_ - start_stats.num_backtrack_operations_;
      delta.num_nonlinear_iters_             = end_stats.num_nonlinear_iters_ - start_stats.num_nonlinear_iters_;
      delta.num_nonlinear_convergence_fails_ = end_stats.num_nonlinear_convergence_fails_
                                               - start_stats.num_nonlinear_convergence_fails_;
      delta.num_nonlinear_step_fails_     = end_stats.num_nonlinear_step_fails_ - start_stats.num_nonlinear_step_fails_;
      delta.num_jacobian_evals_           = end_stats.num_jacobian_evals_ - start_stats.num_jacobian_evals_;
      delta.num_linear_iters_             = end_stats.num_linear_iters_ - start_stats.num_linear_iters_;
      delta.num_linear_convergence_fails_ = end_stats.num_linear_convergence_fails_
                                            - start_stats.num_linear_convergence_fails_;
      delta.num_linear_residual_evals_ = end_stats.num_linear_residual_evals_ - start_stats.num_linear_residual_evals_;
      delta.num_preconditioner_evals_  = end_stats.num_preconditioner_evals_ - start_stats.num_preconditioner_evals_;
      delta.num_preconditioner_solves_ = end_stats.num_preconditioner_solves_ - start_stats.num_preconditioner_solves_;
      delta.num_jtimes_setup_evals_    = end_stats.num_jtimes_setup_evals_ - start_stats.num_jtimes_setup_evals_;
      delta.num_jtimes_evals_          = end_stats.num_jtimes_evals_ - start_stats.num_jtimes_evals_;

      return delta;
    }

    void writeIdaStatsJson(const IdaStatsReport& report, const IdaDiagnosticsOutput& output)
    {
      std::ofstream stream(output.file);
      if (!stream)
      {
        throw std::runtime_error("failed to open IDA stats output file '" + output.file.string() + "'");
      }

      auto json             = idaStatsJson(report.summary);
      json["segment_count"] = report.segments.size();
      json["segments"]      = nlohmann::json::array();
      for (const auto& segment : report.segments)
      {
        json["segments"].push_back(idaStatsSegmentJson(segment));
      }

      const auto& log = output.log.has_value() ? output.log : report.log;
      if (log.has_value())
      {
        json["log"] = {
            {"file", log->file.string()},
            {"level", logLevelName(log->level)}};
      }

      stream << json.dump(2) << '\n';
    }
  } // namespace Sundials
} // namespace AnalysisManager
