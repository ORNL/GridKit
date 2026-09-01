/**
 * @file TwoBusBasic.cpp
 * @author Adam Birchfield (abirchfield@tamu.edu)
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Example running a 2-bus system
 *
 * Simulates a 2-bus system with Genrou 6th order generator model and
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

int main(int argc, const char* argv[])
{
  using namespace GridKit::PhasorDynamics;
  using namespace AnalysisManager::Sundials;

  using scalar_type = double;
  using real_type   = double;
  using index_type  = size_t;

  //
  // Input file
  //

  std::filesystem::path input_file;
  if (argc < 2)
  {
    if (std::filesystem::exists("TwoBusBasic.case.json"))
    {
      input_file = std::filesystem::current_path() / "TwoBusBasic.case.json";
    }
    else
    {
      std::cout << "\n"
                   "ERROR: No input file found or provided.\n"
                   "\n"
                   "Usage:\n"
                   "       TwoBusBasicJson <json-input-file>\n"
                   "\n"
                   "Please provide a JSON input file as a positional command-line \n"
                   "argument.\n"
                   "\n"
                   "By default this example will look for \"TwoBusBasic.json\" in the \n"
                   "current working directory and use that if found.\n"
                   "\n";
      exit(1);
    }
  }
  else
  {
    input_file = argv[1];
  }

  std::cout << "Example: TwoBusBasicJson\n";
  std::cout << "Input file: " << input_file << '\n';

  //
  // Create model data
  //

  auto data = parseSystemModelData(input_file);

  //
  // Instantiate system model
  //

  SystemModel<scalar_type, index_type> sys(data);
  sys.allocate();

  // Get access to the fault and drive it through its status signal
  auto* fault = sys.getBusFault(0);

  scalar_type                         status_value{0.0};
  index_type                          status_index{GridKit::INVALID_INDEX<index_type>};
  SignalNode<scalar_type, index_type> status_node;
  status_node.set(&status_value, &status_index);
  fault->getSignals().attachSignalNode<BusFaultExternalVariables::STATUS>(&status_node);

  // Monitor every quarter cycle at 60 Hz
  real_type dt_monitor = 1.0 / 4.0 / 60.0;

  // Set up simulation
  Ida<scalar_type, size_t> ida(&sys);
  ida.setMaxSteps(10000);
  ida.configureSimulation();

  // Run simulation
  real_type start = static_cast<real_type>(clock());

  // Run for 1s
  ida.initializeSimulation(0.0, false);
  ida.runSimulation(1.0, dt_monitor);

  // Introduce fault and run for the next 0.1s
  status_node.init(1.0);
  ida.initializeSimulation(1.0);
  ida.runSimulation(1.1, dt_monitor);

  // Clear the fault and run until t = 10s.
  status_node.init(0.0);
  ida.initializeSimulation(1.1);
  ida.runSimulation(10.0, dt_monitor);
  real_type stop = static_cast<real_type>(clock());

  std::cout << "\n\nComplete in " << (stop - start) / CLOCKS_PER_SEC << " seconds\n";

  return 0;
}
