#include "PDSim.hpp"

#include <exception>
#include <filesystem>
#include <fstream>
#include <optional>
#include <string>

#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>
#include <GridKit/Solver/Dynamic/IdaDiagnostics.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>

using Log = GridKit::Utilities::Logger;

using namespace GridKit::PhasorDynamics;
using namespace GridKit::Testing;
using namespace AnalysisManager::Sundials;

using scalar_type = double;
using real_type   = double;
using index_type  = size_t;

int main(int argc, const char* argv[])
{
  // Study file
  if (argc < 2)
  {
    Log::error() << "No input file provided" << std::endl;
    std::cout << "\n"
                 "Usage:\n"
                 "       pdsim <json-input-file>\n"
                 "\n"
                 "Please provide a json input file for the study to run.\n"
                 "\n";
    exit(1);
  }

  auto study = parseStudyData(argv[1]);

  // Instantiate system
  SystemModel<scalar_type, index_type> sys(study.model_data);
  sys.allocate();

  real_type dt = study.dt;

  // Set up simulation
  Ida<scalar_type, index_type> ida(&sys);
  IdaStatsRecorder             ida_stats(!study.ida_stats_file.empty());
  ida.configureSimulation();

  // Start timer
  real_type start = static_cast<real_type>(clock());

  using EventType = SystemEvent::Type;

  // Initilize simultation for first run
  ida.initializeSimulation(0.0, false);
  real_type                  curr_time    = 0.0;
  int                        solve_status = 0;
  std::exception_ptr         pending_exception;
  std::optional<std::string> pending_what;
  for (const auto& event : study.events)
  {
    // Run to event time
    int nout = static_cast<int>(std::round((event.time - curr_time) / dt));
    if (nout > 0)
    {
      ida_stats.beginSegment(ida);
      try
      {
        solve_status = ida.runSimulation(event.time, nout);
      }
      catch (const std::exception& ex)
      {
        solve_status      = -1;
        pending_exception = std::current_exception();
        pending_what      = std::string(ex.what());
      }
      // Capture stats even if the segment failed — IdaGetX counters remain valid.
      ida_stats.endSegment(ida, curr_time, event.time, nout);
    }
    if (solve_status != 0)
    {
      break;
    }

    // Set up run for event (to start at event time)
    if (event.type == EventType::FAULT_ON)
    {
      sys.getBusFault(event.element_id)->setStatus(true);
    }
    else if (event.type == EventType::FAULT_OFF)
    {
      sys.getBusFault(event.element_id)->setStatus(false);
    }

    // Re-initialize simulation at event time
    ida.initializeSimulation(event.time, true);
    curr_time = event.time;
  }

  // Run to final time
  int nout = static_cast<int>(std::round((study.tmax - curr_time) / dt));
  if (solve_status == 0 && nout > 0)
  {
    ida_stats.beginSegment(ida);
    try
    {
      solve_status = ida.runSimulation(study.tmax, nout);
    }
    catch (const std::exception& ex)
    {
      solve_status      = -1;
      pending_exception = std::current_exception();
      pending_what      = std::string(ex.what());
    }
    ida_stats.endSegment(ida, curr_time, study.tmax, nout);
  }

  real_type stop = static_cast<real_type>(clock());

  // Stop the variable monitor
  sys.stopMonitor();

  if (!study.ida_stats_file.empty())
  {
    writeIdaStatsJson(ida_stats.report(), {study.ida_stats_file, std::nullopt});
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
  std::string func{"monitor file vs reference file"};
  TestStatus  status{func.c_str()};
  if (!study.output_file.empty() && !study.reference_file.empty())
  {
    auto errorSet = compareCSV(study.output_file, study.reference_file);

    // Print the errors
    errorSet.display();

    // Check against specified tolerance
    status *= errorSet.total.max_error < study.error_tol;

    status.report();
  }

  // Report run time
  std::cout << "\n\nComplete in " << (stop - start) / CLOCKS_PER_SEC << " seconds\n";

  return solve_status == 0 ? status.get() : solve_status;
}
