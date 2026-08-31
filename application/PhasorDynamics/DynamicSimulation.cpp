#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>

#include <GridKit/Model/PhasorDynamics/BusFault/BusFault.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>

#include "AnalysisUtilities.hpp"

using namespace AnalysisManager::Sundials;
using namespace GridKit::PhasorDynamics;
using namespace GridKit::Testing;

using scalar_type = double;
using real_type   = double;
using index_type  = size_t;

int main(int argc, const char* argv[])
{
  using ProfileClock           = std::chrono::steady_clock;
  const auto application_start = ProfileClock::now();

  // Study file
  checkCommandLine(argc, "DynamicSimulation");
  const auto parse_start   = ProfileClock::now();
  auto       study         = parseStudyData(argv[1]);
  const auto parse_seconds = std::chrono::duration<double>(ProfileClock::now() - parse_start).count();

  // Instantiate system
  const auto                           model_construction_start = ProfileClock::now();
  SystemModel<scalar_type, index_type> sys(study.model_data);
  const auto                           model_construction_seconds = std::chrono::duration<double>(ProfileClock::now() - model_construction_start).count();

  const auto model_allocation_start = ProfileClock::now();
  sys.allocate();
  const auto model_allocation_seconds = std::chrono::duration<double>(ProfileClock::now() - model_allocation_start).count();

  // Set up simulation
  const auto                   ida_configuration_start = ProfileClock::now();
  Ida<scalar_type, index_type> ida(&sys);
  ida.setTolerance(study.rel_tol, study.abs_tol);
  ida.setFixedStep(study.dt_fixed);
  ida.setMaxSteps(study.max_steps);
  ida.setConsistentICType(study.consistent_ic_type);
  ida.configureSimulation();
  const auto ida_configuration_seconds = std::chrono::duration<double>(ProfileClock::now() - ida_configuration_start).count();

  const auto differential_variables = static_cast<index_type>(std::count(sys.tag().begin(), sys.tag().end(), true));
  std::cout << "\nGRIDKIT_SYSTEM_BEGIN\n"
            << "buses=" << study.model_data.bus.size() << '\n'
            << "branches=" << study.model_data.branch.size() << '\n'
            << "regca=" << study.model_data.regca.size() << '\n'
            << "loadz=" << study.model_data.loadz.size() << '\n'
            << "loadzip=" << study.model_data.loadzip.size() << '\n'
            << "genrou=" << study.model_data.genrou.size() << '\n'
            << "gensal=" << study.model_data.gensal.size() << '\n'
            << "genclassical=" << study.model_data.genclassical.size() << '\n'
            << "states=" << sys.size() << '\n'
            << "differential_variables=" << differential_variables << '\n'
            << "algebraic_variables=" << sys.size() - differential_variables << '\n'
            << "jacobian_nnz=" << sys.nnz() << '\n'
            << "GRIDKIT_SYSTEM_END\n";

  // Allocation performs an untimed sparsity-discovery evaluation. Exclude it
  // from the simulation profile so every metric has the same timed boundary.
  sys.resetPerformanceStats();
  ida.resetPerformanceStats();

  // Start timer
  const auto wall_start = ProfileClock::now();
  real_type  start      = static_cast<real_type>(clock());

  using EventType = SystemEvent::Type;

  // Initilize simultation for first run
  auto      dt_monitor = study.dt_monitor;
  real_type final_time = study.tmax;
  IdaStats  stats;
  ida.initializeSimulation(0.0);
  for (const auto& event : study.events)
  {
    // Run to event time
    ida.runSimulation(event.time, dt_monitor);
    stats += ida.getStats();

    // Set up run for event (to start at event time)
    switch (event.type)
    {
    case EventType::FAULT_ON:
      sys.getBusFault(event.element_id)->setStatus(true);
      break;
    case EventType::FAULT_OFF:
      sys.getBusFault(event.element_id)->setStatus(false);
      break;
    }

    // Re-initialize simulation at event time
    ida.initializeSimulation(event.time);
  }

  // Run to final time
  ida.runSimulation(final_time, dt_monitor);
  stats += ida.getStats();

  real_type  stop         = static_cast<real_type>(clock());
  const auto wall_seconds = std::chrono::duration<double>(ProfileClock::now() - wall_start).count();

  // Stop the variable monitor
  const auto postprocess_start = ProfileClock::now();
  sys.stopMonitor();

  // Generate aggregate errors comparing variable output to reference solution
  TestStatus status              = checkErrors(study);
  const auto postprocess_seconds = std::chrono::duration<double>(ProfileClock::now() - postprocess_start).count();
  const auto application_seconds = std::chrono::duration<double>(ProfileClock::now() - application_start).count();

  // Report run time
  std::cout << "\n\nComplete in " << (stop - start) / CLOCKS_PER_SEC << " seconds\n";
  std::cout << std::fixed << std::setprecision(9)
            << "\nAPPLICATION_PROFILE_BEGIN\n"
            << "study_parse_seconds=" << parse_seconds << '\n'
            << "model_construction_seconds=" << model_construction_seconds << '\n'
            << "model_allocation_seconds=" << model_allocation_seconds << '\n'
            << "ida_configuration_seconds=" << ida_configuration_seconds << '\n'
            << "simulation_wall_seconds=" << wall_seconds << '\n'
            << "postprocess_seconds=" << postprocess_seconds << '\n'
            << "application_wall_seconds=" << application_seconds << '\n'
            << "APPLICATION_PROFILE_END\n"
            << '\n'
            << stats.report() << '\n';
  ida.printPerformanceStats();
  sys.printPerformanceStats();

  return status.get();
}
