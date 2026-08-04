#include <exception>
#include <filesystem>
#include <fstream>
#include <optional>
#include <string>

#include <GridKit/Model/PhasorDynamics/Branch/Branch.hpp>
#include <GridKit/Model/PhasorDynamics/BusBase.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>
#include <GridKit/Solver/Dynamic/IdaDiagnostics.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>

#include "AnalysisUtilities.hpp"

using Log = GridKit::Utilities::Logger;

using namespace GridKit::PhasorDynamics;
using namespace GridKit::Testing;
using namespace AnalysisManager::Sundials;

using scalar_type = double;
using real_type   = double;
using index_type  = size_t;

int runApplication(int argc, const char* argv[])
{
  // Print summaries, such as the run time, without lowering a higher verbosity
  Log::raiseVerbosity(Log::SUMMARY);

  // Study file
  checkCommandLine(argc, "DynamicSimulation");
  auto study = parseStudyData(argv[1]);

  using EventType = SystemEvent::Type;

  for (const auto& event : study.events)
  {
    if ((event.type == EventType::FAULT_ON || event.type == EventType::FAULT_OFF)
        && event.bus == GridKit::INVALID_INDEX<index_type>)
    {
      Log::error() << "Fault event at time " << event.time
                   << " requires a bus number\n";
      return 1;
    }
  }

  // Instantiate system
  SystemModel<scalar_type, index_type> sys(study.model_data);
  sys.allocate();

  // Set up simulation
  Ida<scalar_type, index_type> ida(&sys);
  IdaStatsRecorder             ida_stats_recorder(!study.ida_stats.empty());
  IdaStepHistoryRecorder       ida_steps_recorder(!study.ida_steps.empty());
  ida.setTolerance(study.rel_tol, study.abs_tol);
  ida.setFixedStep(study.dt_fixed);
  ida.setMaxSteps(study.max_steps);
  ida.setMaxOrder(study.max_order);
  ida.setConsistentICType(study.consistent_ic_type);
  ida.configureSimulation();
  if (study.ida_max_dt.has_value())
  {
    if (study.ida_max_dt.value() <= 0.0)
    {
      Log::error() << "ida_max_dt must be positive" << std::endl;
      return 1;
    }
    ida.setMaxStep(static_cast<real_type>(study.ida_max_dt.value()));
  }

  // Start timer
  real_type start = static_cast<real_type>(clock());

  // Initialize simulation for first run
  // A negative dt_monitor in the study file requests output at
  // solver-selected steps, passed to IDA as an empty monitor interval.
  const std::optional<real_type> dt_monitor = study.dt_monitor < 0.0
                                                  ? std::optional<real_type>{}
                                                  : std::optional<real_type>{study.dt_monitor};
  real_type                      final_time = study.tmax;
  ida.initializeSimulation(0.0);
  real_type                  curr_time    = 0.0;
  int                        solve_status = 0;
  std::exception_ptr         pending_exception;
  std::optional<std::string> pending_what;
  auto                       run_segment = [&](real_type start_time, real_type end_time)
  {
    int                                                         output_steps = 0;
    std::optional<Ida<scalar_type, index_type>::OutputCallback> output_callback;
    if (ida_stats_recorder.enabled() || ida_steps_recorder.enabled())
    {
      output_callback = [&](real_type)
      { ++output_steps; };
    }
    ida_stats_recorder.beginSegment(ida);
    ida_steps_recorder.beginSegment(ida, start_time, end_time);
    try
    {
      if (ida_steps_recorder.enabled())
      {
        solve_status = ida.runSimulationWithStepHistory(end_time, dt_monitor, [&](const IdaStats& stats)
                                                        { ida_steps_recorder.recordStep(stats); },
                                                        output_callback);
      }
      else
      {
        solve_status = ida.runSimulation(end_time, dt_monitor, output_callback);
      }
    }
    catch (const std::exception& ex)
    {
      solve_status      = -1;
      pending_exception = std::current_exception();
      pending_what      = std::string(ex.what());
    }
    // Capture diagnostics even if the segment failed — IdaGetX counters remain valid.
    ida_steps_recorder.endSegment(ida, output_steps);
    ida_stats_recorder.endSegment(ida, start_time, end_time, output_steps);
  };

  for (const auto& event : study.events)
  {
    // Run to event time
    if (event.time > curr_time)
    {
      run_segment(curr_time, event.time);
    }
    if (solve_status != 0)
    {
      break;
    }

    // Set up run for event (to start at event time)
    int ret = 0;

    switch (event.type)
    {
    case EventType::FAULT_ON:
      ret = sys.getBus(event.bus)->setFault(true, event.R, event.X);
      break;
    case EventType::FAULT_OFF:
      ret = sys.getBus(event.bus)->clearFault();
      break;
    case EventType::IN_SERVICE:
      sys.getBranch(event.device)->setInService(true);
      break;
    case EventType::OUT_OF_SERVICE:
      sys.getBranch(event.device)->setInService(false);
      break;
    }

    if (ret != 0)
    {
      Log::error() << "Fault event failed for bus " << event.bus
                   << " at time " << event.time << '\n';
      return 1;
    }

    // Re-initialize simulation at event time
    ida.initializeSimulation(event.time);
    curr_time = event.time;
  }

  // Run to final time
  if (solve_status == 0 && final_time > curr_time)
  {
    run_segment(curr_time, final_time);
  }

  real_type stop = static_cast<real_type>(clock());

  // Stop the variable monitor
  sys.stopMonitor();

  if (!study.ida_stats.empty())
  {
    writeIdaStatsJson(ida_stats_recorder.report(), {study.ida_stats, std::nullopt});
  }
  if (!study.ida_steps.empty())
  {
    writeIdaStepHistoryJson(ida_steps_recorder.report(), {study.ida_steps, std::nullopt});
  }

  // Preserve original behaviour: surface the SUNDIALS exception (and its message)
  // by re-throwing now that diagnostics have been flushed to disk.
  if (pending_exception)
  {
    if (pending_what.has_value())
    {
      Log::error() << *pending_what << std::endl;
    }
    std::rethrow_exception(pending_exception);
  }

  // Generate aggregate errors comparing variable output to reference solution
  TestStatus status = checkErrors(study);

  // Report run time
  Log::summary() << "Complete in " << (stop - start) / CLOCKS_PER_SEC << " seconds\n";

  return solve_status == 0 ? status.get() : solve_status;
}

int main(int argc, const char* argv[])
{
  try
  {
    return runApplication(argc, argv);
  }
  catch (const std::exception& error)
  {
    Log::error() << "DynamicSimulation failed: " << error.what() << '\n';
  }

  return 1;
}
