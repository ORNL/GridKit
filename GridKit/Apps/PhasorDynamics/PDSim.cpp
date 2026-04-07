#include "PDSim.hpp"

#include <filesystem>
#include <fstream>

#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>
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
  ida.configureSimulation();

  // Run simulation, output each `dt` interval
  real_type start = static_cast<real_type>(clock());

  ida.initializeSimulation(0.0, false);

  real_type curr_time = 0.0;
  for (const auto& cue : study.schedule)
  {
    // Run to scheduled time
    int nout = static_cast<int>(std::round((cue.time - curr_time) / dt));
    ida.runSimulation(cue.time, nout);

    // Execute action
    const auto& ev = study.event_map.at(cue.event);
    if (ev.type == "bus_fault")
    {
      sys.getBusFault(study.fault_map.at(cue.event))->setStatus(cue.action == "on");
    }

    ida.initializeSimulation(cue.time, false);
    curr_time = cue.time;
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

    status *= errorSet.total.max_error < study.error_tol;

    status.report();
  }

  std::cout << "\n\nComplete in " << (stop - start) / CLOCKS_PER_SEC << " seconds\n";

  return status.get();
}
