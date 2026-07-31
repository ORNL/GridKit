#include <filesystem>
#include <fstream>

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

  // Set up simulation
  Ida<scalar_type, index_type> ida(&sys);
  ida.setTolerance(study.rel_tol, study.abs_tol);
  ida.setFixedStep(study.dt_fixed);
  ida.configureSimulation();

  // Start timer
  real_type start = static_cast<real_type>(clock());

  // Initilize simultation for first run
  auto      dt_monitor = study.dt_monitor;
  real_type final_time = study.tmax;
  ida.initializeSimulation(0.0);
  for (std::size_t i = 0; i < study.events.size();)
  {
    const real_type event_time = study.events[i].time;

    // Run to event time
    ida.runSimulation(event_time, dt_monitor);

    // Apply every input change scheduled for this time before reinitializing.
    do
    {
      const auto& event = study.events[i];
      sys.setInput(event.device_id, event.input, event.value);
      ++i;
    } while (i < study.events.size()
             && study.events[i].time == event_time);

    // Re-initialize simulation at event time
    ida.initializeSimulation(event_time);
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
