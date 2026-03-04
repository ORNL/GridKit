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
#include <iostream>

#include <GridKit/Model/PhasorDynamics/ComponentLibrary.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModelData.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/CliArgs/CliArgs.hpp>

#define ERROR_TOL 1.0e-4

using scalar_type = double;
using real_type   = double;
using index_type  = size_t;

int main(int argc, const char* argv[])
{
  using namespace GridKit::PhasorDynamics;
  using namespace AnalysisManager::Sundials;
  using namespace GridKit::Utilities;
  using namespace GridKit::Testing;

  CliArgs args{
      {
          .name     = {"--case", "-c"},
          .help     = "JSON file describing the system configuration",
          .required = false,
          .defaults = "ThreeBusBasic.case.json",
      },

      {
          .name     = {"--compare", "-r"},
          .help     = "Two CSV files to compare:\n"
                      "<expected-output-file> <reference-file>",
          .nargs    = 2,
          .defaults = {"mon.csv", "ThreeBusBasic.ref.csv"},
      },

      {
          .name     = {"--tolerance", "-t"},
          .help     = "Allowable maximum error between compared files",
          .type     = ArgType::Real,
          .defaults = 1.0e-4,
      }};

  args.parseArgs(argc, argv);

  // Input file
  auto input_file = args["case"]();
  std::cout << "Example: ThreeBusBasicJson\n";
  std::cout << "Input file: " << input_file << '\n';

  // Create model data
  auto data = parseSystemModelData(input_file);

  // Instantiate system
  SystemModel<scalar_type, index_type> sys(data);
  sys.allocate();

  // Get access to fault 0
  auto* fault = sys.getBusFault(0);

  real_type dt = 1.0 / 4.0 / 60.0;

  // Set up simulation
  Ida<scalar_type, index_type> ida(&sys);
  ida.configureSimulation();

  // Run simulation, output each `dt` interval
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

  // Stop the variable monitor
  sys.stopMonitor();

  // Generate aggregate errors comparing variable output to reference solution
  TestStatus status{__func__};
  if (!args["compare"].empty())
  {
    const auto& [out_file, ref_file] = args["compare"].as<2>();

    auto errorSet = GridKit::Testing::compareCSV(out_file, ref_file);

    // Print the errors
    errorSet.display();

    auto error_allowed = args["tolerance"].as<real_type>();

    status *= errorSet.total.max_error < error_allowed;

    status.report();
  }

  std::cout << "\n\nComplete in " << (stop - start) / CLOCKS_PER_SEC << " seconds\n";

  return status.get();
}
