/**
 * @file texas.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @brief Example running the ACTIVSg2000 Texas Synthetic Case
 *
 * Simulates a 2000-bus Texas case with Genrou 6th order generator model.
 *
 */

#include <cmath>
#include <ctime>
#include <exception>
#include <filesystem>
#include <iostream>
#include <string>

#include <GridKit/Model/PhasorDynamics/ComponentLibrary.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>

int main(int argc, const char* argv[])
{
  using namespace GridKit::PhasorDynamics;
  using namespace AnalysisManager::Sundials;

  using scalar_type = double;
  using real_type   = double;
  using index_type  = size_t;

  // Read input JSON file.
  std::filesystem::path input_file;
  if (argc < 2)
  {
    if (std::filesystem::exists("texas.json"))
    {
      input_file = std::filesystem::current_path() / "texas.json";
    }
    else
    {
      std::cout << "\n"
                   "ERROR: No input file found or provided.\n"
                   "\n"
                   "Usage:\n"
                   "       texas <json-input-file>\n"
                   "\n"
                   "Please provide a JSON input file as a positional command-line \n"
                   "argument.\n"
                   "\n"
                   "By default this example will look for \"texas.json\" in the \n"
                   "current working directory and use that if found.\n"
                   "\n";
      return 1;
    }
  }
  else
  {
    input_file = argv[1];
  }

  std::cout << "Example: texas\n";
  std::cout << "Input file: " << input_file << '\n';

  // Parse model data.
  auto data = parseSystemModelData(input_file);

  // Instantiate system model.
  SystemModel<scalar_type, index_type> sys(data);
  sys.allocate();

  // Set time step to 1/4 of a 60Hz cycle.
  real_type dt = 1.0 / 4.0 / 60.0;

  // Set up simulation.
  Ida<scalar_type, index_type> ida(&sys);
  ida.configureSimulation();

  real_type   start = static_cast<real_type>(clock());
  std::string phase = "setup";

  // Get access to fault 0.
  auto* fault = sys.getBusFault(0);
  if (fault == nullptr)
  {
    std::cerr << "ERROR: Texas input is missing BusFault component at index 0.\n";
    sys.stopMonitor();
    return 2;
  }

  try
  {
    phase = "pre-fault";
    std::cout << "[pre-fault] t = 0.0 -> 10.0\n";
    ida.initializeSimulation(0.0, false);
    int nout = static_cast<int>(std::round((10.0 - 0.0) / dt));
    ida.runSimulation(10.0, nout);

    phase = "fault";
    std::cout << "[fault] t = 10.0 -> 10.15\n";
    fault->setStatus(true);
    ida.initializeSimulation(10.0, false);
    nout = static_cast<int>(std::round((10.15 - 10.0) / dt));
    ida.runSimulation(10.15, nout);

    phase = "post-fault";
    std::cout << "[post-fault] t = 10.15 -> 20.0\n";
    fault->setStatus(false);
    ida.initializeSimulation(10.15, false);
    nout = static_cast<int>(std::round((20.0 - 10.15) / dt));
    ida.runSimulation(20.0, nout);
  }
  catch (const std::exception& ex)
  {
    std::cerr << "ERROR: Texas simulation failed during '" << phase << "' phase.\n"
              << "Reason: " << ex.what() << '\n';
    sys.stopMonitor();
    return 2;
  }
  catch (...)
  {
    std::cerr << "ERROR: Texas simulation failed during '" << phase
              << "' phase with unknown exception.\n";
    sys.stopMonitor();
    return 2;
  }

  real_type stop = static_cast<real_type>(clock());

  // Flush monitor output.
  sys.stopMonitor();

  std::cout << "\n\nComplete in " << (stop - start) / CLOCKS_PER_SEC << " seconds\n";

  return 0;
}
