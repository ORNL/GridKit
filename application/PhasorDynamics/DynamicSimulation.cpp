#include <filesystem>
#include <fstream>
#include <vector>

#include <GridKit/Model/PhasorDynamics/BusFault/BusFault.hpp>
#include <GridKit/Model/PhasorDynamics/SignalNode/SignalNode.hpp>
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

  // Instantiate system
  SystemModel<scalar_type, index_type> sys(study.model_data);
  sys.allocate();

  // Status signals driving the faults
  const auto                                       n_faults = study.model_data.bus_fault.size();
  std::vector<scalar_type>                         status_values(n_faults, 0.0);
  std::vector<index_type>                          status_indices(n_faults, GridKit::INVALID_INDEX<index_type>);
  std::vector<SignalNode<scalar_type, index_type>> status_nodes(n_faults);
  for (std::size_t i = 0; i < n_faults; ++i)
  {
    status_nodes[i].set(&status_values[i], &status_indices[i]);
    sys.getBusFault(i)->getSignals().attachSignalNode<BusFaultExternalVariables::STATUS>(&status_nodes[i]);
  }

  // Set up simulation
  Ida<scalar_type, index_type> ida(&sys);
  ida.setTolerance(study.rel_tol, study.abs_tol);
  ida.setFixedStep(study.dt_fixed);
  ida.setMaxSteps(study.max_steps);
  ida.setConsistentICType(study.consistent_ic_type);
  ida.configureSimulation();

  // Start timer
  real_type start = static_cast<real_type>(clock());

  using EventType = SystemEvent::Type;

  // Initilize simultation for first run
  auto      dt_monitor = study.dt_monitor;
  real_type final_time = study.tmax;
  ida.initializeSimulation(0.0);
  for (const auto& event : study.events)
  {
    // Run to event time
    ida.runSimulation(event.time, dt_monitor);

    // Set up run for event (to start at event time)
    switch (event.type)
    {
    case EventType::FAULT_ON:
      status_nodes[event.element_id].init(1.0);
      break;
    case EventType::FAULT_OFF:
      status_nodes[event.element_id].init(0.0);
      break;
    }

    // Re-initialize simulation at event time
    ida.initializeSimulation(event.time);
  }

  // Run to final time
  ida.runSimulation(final_time, dt_monitor);

  real_type stop = static_cast<real_type>(clock());

  // Stop the variable monitor
  sys.stopMonitor();

  // Generate aggregate errors comparing variable output to reference solution
  TestStatus status = checkErrors(study);

  // Report run time
  std::cout << "\n\nComplete in " << (stop - start) / CLOCKS_PER_SEC << " seconds\n";

  return status.get();
}
