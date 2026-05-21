#include <filesystem>
#include <fstream>

#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>

#include "AnalysisUtilities.hpp"

using namespace GridKit::PhasorDynamics;
using namespace GridKit::Testing;
using namespace AnalysisManager::Sundials;

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

  real_type dt = study.dt;

  // Set up simulation
  Ida<scalar_type, index_type> ida(&sys);
  ida.configureSimulation();

  // Start timer
  real_type start = static_cast<real_type>(clock());

  using EventType = SystemEvent::Type;

  // Initilize simultation for first run
  ida.initializeSimulation(0.0, false);
  real_type curr_time = 0.0;
  for (const auto& event : study.events)
  {
    // Run to event time
    int nout = static_cast<int>(std::round((event.time - curr_time) / dt));
    ida.runSimulation(event.time, nout);

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
  ida.runSimulation(study.tmax, nout);

  real_type stop = static_cast<real_type>(clock());

  // Stop the variable monitor
  sys.stopMonitor();

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

  return status.get();
}
