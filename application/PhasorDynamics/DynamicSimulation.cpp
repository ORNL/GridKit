#include <filesystem>
#include <fstream>
#include <stdexcept>

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
  // Study file
  checkCommandLine(argc, "DynamicSimulation");
  auto study = parseStudyData(argv[1]);

  if (study.fault_bus)
  {
    if (study.model_data.bus_fault.empty())
    {
      throw std::invalid_argument("fault_bus requires a BusFault in the system model");
    }

    bool bus_found = false;
    for (const auto& bus : study.model_data.bus)
    {
      if (bus.bus_id == *study.fault_bus)
      {
        bus_found = true;
        break;
      }
    }
    if (!bus_found)
    {
      throw std::invalid_argument("fault_bus does not identify a bus in the system model");
    }

    study.model_data.bus_fault.front().buses[BusFaultBuses::bus] = *study.fault_bus;
  }

  // Instantiate system
  SystemModel<scalar_type, index_type> sys(study.model_data);
  sys.allocate();

  std::cout << "\nGRIDKIT_SYSTEM_BEGIN\n"
            << "buses=" << study.model_data.bus.size() << '\n'
            << "states=" << sys.size() << '\n'
            << "jacobian_nnz=" << sys.nnz() << '\n'
            << "GRIDKIT_SYSTEM_END\n";

  // Set up simulation
  Ida<scalar_type, index_type> ida(&sys);
  ida.setOptions(study.ida);
  ida.setConsistentICType(study.consistent_ic_type);
  ida.configureSimulation();
  ida.enableStepTrace(!study.step_trace_file.empty());

  // Start timer
  real_type start = static_cast<real_type>(clock());

  using EventType = SystemEvent::Type;

  // Initilize simultation for first run
  auto      dt_monitor = study.dt_monitor;
  real_type final_time = study.tmax;
  IdaStats  stats;
  int       segment = 0;
  ida.initializeSimulation(0.0);
  for (const auto& event : study.events)
  {
    // Run to event time
    ida.setTraceSegment(segment++);
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
  ida.setTraceSegment(segment);
  ida.runSimulation(final_time, dt_monitor);
  stats += ida.getStats();

  real_type stop = static_cast<real_type>(clock());

  // Stop the variable monitor
  sys.stopMonitor();

  if (!study.step_trace_file.empty())
  {
    writeStepTrace(study.step_trace_file, ida.getStepTrace());
  }

  // Generate aggregate errors comparing variable output to reference solution
  TestStatus status = checkErrors(study);

  // Report run time
  std::cout << "\n\nComplete in " << (stop - start) / CLOCKS_PER_SEC << " seconds\n";
  std::cout << '\n'
            << stats.report() << '\n';
  ida.printPerformanceStats();
  sys.printResidualPerformanceStats();

  return status.get();
}
