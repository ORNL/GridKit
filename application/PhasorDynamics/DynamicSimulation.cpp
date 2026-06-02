#include <cmath>
#include <ctime>
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

TestStatus runEventStudy(StudyData& study)
{
  // Instantiate system
  SystemModel<scalar_type, index_type> sys(study.model_data);
  sys.allocate();

  // Set up simulation
  Ida<scalar_type, index_type> ida(&sys);
  ida.configureSimulation();

  // Start timer
  real_type start = static_cast<real_type>(clock());

  using EventType = SystemEvent::Type;

  // Initilize simultation for first run
  real_type dt         = study.dt;
  real_type final_time = study.tmax;
  real_type curr_time  = 0.0;
  ida.initializeSimulation(0.0, false);
  for (const auto& event : study.events)
  {
    // Run to event time
    int nout = static_cast<int>(std::round((event.time - curr_time) / dt));
    ida.runSimulation(event.time, nout);

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
    ida.initializeSimulation(event.time, true);
    curr_time = event.time;
  }

  // Run to final time
  int nout = static_cast<int>(std::round((final_time - curr_time) / dt));
  ida.runSimulation(final_time, nout);

  real_type stop = static_cast<real_type>(clock());

  // Stop the variable monitor
  sys.stopMonitor();

  // Generate aggregate errors comparing variable output to reference solution
  TestStatus status = checkErrors(study);

  // Report run time
  std::cout << "\n\nComplete in " << (stop - start) / CLOCKS_PER_SEC << " seconds\n";

  return status;
}

TestStatus runForcedOscillationStudy(StudyData& study)
{
  auto model_data = applyInjections(study.base_model_data, study.forced_oscillations);

  // Instantiate system
  SystemModel<scalar_type, index_type> sys(model_data);
  sys.allocate();

  // Set up simulation
  Ida<scalar_type, index_type> ida(&sys);
  ida.configureSimulation();
  if (study.max_step.has_value())
  {
    ida.setMaxStep(*study.max_step);
  }

  // Start timer
  real_type start = static_cast<real_type>(clock());

  // Initialize simulation once and run directly to final time
  real_type final_time = study.tmax;
  ida.initializeSimulation(0.0, false);
  int nout = static_cast<int>(std::round(final_time / study.dt));
  ida.runSimulation(final_time, nout);

  real_type stop = static_cast<real_type>(clock());

  // Stop the variable monitor
  sys.stopMonitor();

  // Generate aggregate errors comparing variable output to reference solution
  TestStatus status = checkErrors(study);

  // Report run time
  std::cout << "\n\nComplete in " << (stop - start) / CLOCKS_PER_SEC << " seconds\n";

  return status;
}

int main(int argc, const char* argv[])
{
  // Study file
  checkCommandLine(argc, "DynamicSimulation");
  auto study = parseStudyData(argv[1]);

  TestStatus status = true;
  if (study.isForcedOscillationStudy())
  {
    status = runForcedOscillationStudy(study);
  }
  else
  {
    status = runEventStudy(study);
  }

  return status.get();
}
