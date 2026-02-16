/**
 * @file ThreeBusBasicJson.cpp
 * @author Adam Birchfield (abirchfield@tamu.edu)
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Example running a 3-bus system
 *
 * Simulates a 3-bus system with two Genrou 6th order generator models and
 * compares results with data generated for the same system by Poweworld.
 *
 */
#include <cmath>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <vector>

#include <GridKit/Model/PhasorDynamics/ComponentLibrary.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>
#include <GridKit/Utilities/CliOptions/CliOptions.hpp>
#include <GridKit/Utilities/TestHelpers.hpp>
#include <GridKit/Utilities/Testing.hpp>

#define ERROR_TOL 1.0e-4

using scalar_type = double;
using real_type   = double;
using index_type  = size_t;

// NOTES:
// Write function to compare to CSV files
//  compare number of lines and number of items per line

int main(int argc, const char* argv[])
{
  using namespace GridKit::PhasorDynamics;
  using namespace AnalysisManager::Sundials;
  using namespace GridKit::Utilities;

  auto error_allowed = static_cast<real_type>(1e-4);

  //
  // Input file
  //

  std::filesystem::path input_file;
  if (argc < 2)
  {
    if (std::filesystem::exists("ThreeBusBasic.json"))
    {
      input_file = std::filesystem::current_path() / "ThreeBusBasic.json";
    }
    else
    {
      std::cout << "\n"
                   "ERROR: No input file found or provided.\n"
                   "\n"
                   "Usage:\n"
                   "       ThreeBusBasicJson <json-input-file>\n"
                   "\n"
                   "Please provide a JSON input file as a positional command-line \n"
                   "argument.\n"
                   "\n"
                   "By default this example will look for \"ThreeBusBasic.json\" in the \n"
                   "current working directory and use that if found.\n"
                   "\n";
      exit(1);
    }
  }
  else
  {
    input_file = argv[1];
  }

  std::cout << "Example: ThreeBusBasicJson\n";
  std::cout << "Input file: " << input_file << '\n';

  //
  // Create model data
  //

  auto data = parseSystemModelData(input_file);

  //
  // Instantiate system
  //

  SystemModel<scalar_type, index_type> sys(data);
  sys.allocate();

  // Get access to fault 0
  auto* fault = sys.getBusFault(0);

  real_type dt = 1.0 / 4.0 / 60.0;

  // Set up simulation
  Ida<scalar_type, index_type> ida(&sys);
  ida.configureSimulation();

  // Run simulation, output each `dt` interval
  scalar_type start = static_cast<scalar_type>(clock());
  ida.initializeSimulation(0.0, false);

  // Run for 1s
  int nout = static_cast<int>(std::round((1.0 - 0.0) / dt));
  ida.runSimulation(1.0, nout);

  // Introduce fault to ground and run for 0.1s
  fault->setStatus(true);
  ida.initializeSimulation(1.0, false);
  nout = static_cast<int>(std::round((1.1 - 1.0) / dt));
  ida.runSimulation(1.1, nout);

  // Clear fault and run until t = 10s.
  fault->setStatus(false);
  ida.initializeSimulation(1.1, false);
  nout = static_cast<int>(std::round((10.0 - 1.1) / dt));
  ida.runSimulation(10.0, nout);
  double stop = static_cast<double>(clock());

  // Stop the variable monitor
  sys.stopMonitor();

  // Generate aggregate errors comparing variable output to reference solution
  auto errorSet =
      GridKit::Testing::compareCSV("mon.csv", "ThreeBusBasic.ref.csv");

  // Print the errors
  errorSet.display();

  std::cout << "\n\nComplete in " << (stop - start) / CLOCKS_PER_SEC << " seconds\n";

  return errorSet.total.max_error < error_allowed ? 0 : 1;
}
