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
                 "       pdsim <solver-file.solver.json>\n"
                 "\n"
                 "Please provide a JSON solver file for the study to run.\n"
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

  // Start timer
  real_type start = static_cast<real_type>(clock());

  // Initilize simultation for first run
  ida.initializeSimulation(0.0, false);
  real_type curr_time = 0.0;

  size_t i = 0;
  while (i < study.schedule.size())
  {
    real_type t_cue = study.schedule[i].time;

    // Run to cue time
    int nout = static_cast<int>(std::round((t_cue - curr_time) / dt));
    ida.runSimulation(t_cue, nout);

    // Batch simultaneous cues into one re-init so IDA never sees a
    // zero-length segment between back-to-back applies at the same time.
    while (i < study.schedule.size() && study.schedule[i].time == t_cue)
    {
      sys.cue(study.schedule[i].target, study.schedule[i].action);
      ++i;
    }

    // Re-initialize simulation at cue time
    ida.initializeSimulation(t_cue, true);
    curr_time = t_cue;
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
