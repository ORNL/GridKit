/**
 * @file hawaiiJson.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @author Adam Birchfield (abirchfield@tamu.edu)
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Example running the 37-Bus Hawaii Synthetic Case
 *
 * Simulates a Hawaii Case with Genrou 6th order generator model and
 * compares results with data generated for the same system by Poweworld.
 *
 */
#include <ctime>
#include <filesystem>
#include <iostream>

#include <GridKit/Model/PhasorDynamics/ComponentLibrary.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>
#include <GridKit/Testing/Testing.hpp>

#include "hawaii.hpp"

int main(int argc, const char* argv[])
{
  using namespace GridKit::PhasorDynamics;
  using namespace AnalysisManager::Sundials;

  using scalar_type = double;
  using real_type   = double;
  using index_type  = size_t;

  // Read Input JSON File
  std::filesystem::path input_file;
  if (argc < 2)
  {
    if (std::filesystem::exists("hawaii.json"))
    {
      input_file = std::filesystem::current_path() / "hawaii.json";
    }
    else
    {
      std::cout << "\n"
                   "ERROR: No input file found or provided.\n"
                   "\n"
                   "Usage:\n"
                   "       hawaiiJson <json-input-file>\n"
                   "\n"
                   "Please provide a JSON input file as a positional command-line \n"
                   "argument.\n"
                   "\n"
                   "By default this example will look for \"hawaii.json\" in the \n"
                   "current working directory and use that if found.\n"
                   "\n";
      exit(1);
    }
  }
  else
  {
    input_file = argv[1];
  }

  std::cout << "Example: hawaiiJson\n";
  std::cout << "Input file: " << input_file << '\n';

  // Parse Model Data
  auto data = parseSystemModelData(input_file);

  //  Instantiate System Model
  SystemModel<scalar_type, index_type> sys(data);
  sys.allocate();

  // Get access to fault 0
  auto* fault = sys.getBusFault(0);

  // NOTE Now we try and run the case.
  // Fails for now but left here for future

  // Set time step to 1/4 of a 60Hz cycle
  real_type dt = 1.0 / 4.0 / 60.0;

  // Set up simulation
  Ida<scalar_type, index_type> ida(&sys);
  ida.configureSimulation();

  // Initialize simulation
  real_type start = static_cast<real_type>(clock());
  ida.initializeSimulation(0.0, true);

  // Run for 1s
  int nout = static_cast<int>(std::round((1.0 - 0.0) / dt));
  ida.runSimulation(1.0, nout);

  // Introduce fault to ground and run for 0.1s
  fault->setStatus(true);
  ida.initializeSimulation(1.0, true);
  nout = static_cast<int>(std::round((1.1 - 1.0) / dt));
  ida.runSimulation(1.1, nout);

  // Clear fault and run until t = 10s.
  fault->setStatus(false);
  ida.initializeSimulation(1.1, true);
  nout = static_cast<int>(std::round((10.0 - 1.1) / dt));
  ida.runSimulation(10.0, nout);
  real_type stop = static_cast<real_type>(clock());

  std::cout << "\n\nComplete in " << (stop - start) / CLOCKS_PER_SEC << " seconds\n";

  int status = 0;
  return status;
}
