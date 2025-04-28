/**
 * @file example1.cpp
 * @author Adam Birchfield (abirchfield@tamu.edu)
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Example running a 2-bus system
 *
 * Simulates a 2-bus system with Genrou 6th order generator model and
 * compares results with data generated for the same system by Poweworld.
 *
 */
#include "example1.hpp"

#include <ctime>
#include <iostream>

#include <Model/PhasorDynamics/Branch/Branch.hpp>
#include <Model/PhasorDynamics/Bus/Bus.hpp>
#include <Model/PhasorDynamics/Bus/BusInfinite.hpp>
#include <Model/PhasorDynamics/BusFault/BusFault.hpp>
#include <Model/PhasorDynamics/SynchronousMachine/GENROUwS/Genrou.hpp>
#include <Model/PhasorDynamics/SystemModel.hpp>
#include <Solver/Dynamic/Ida.hpp>
#include <Utilities/Testing.hpp>

int main()
{
  using namespace GridKit::PhasorDynamics;
  using namespace AnalysisManager::Sundials;

  std::cout << "Example 1 version 2\n";

  /* Create model parts */
  SystemModel<double, size_t> sys;
  Bus<double, size_t>         bus1(0.9949877346411762, 0.09999703952427966);
  BusInfinite<double, size_t> bus2(1.0, 0.0);
  Branch<double, size_t>      branch(&bus1, &bus2, 0, 0.1, 0, 0);
  BusFault<double, size_t>    fault(&bus1, 0, 1e-3, 0);

  Genrou<double, size_t> gen(&bus1,
                             1,
                             1.,
                             0.05013,
                             3.,
                             0.,
                             0.,
                             7.,
                             .04,
                             .05,
                             .75,
                             2.1,
                             0.2,
                             0.18,
                             0.5,
                             0.5,
                             0.18,
                             0.15,
                             0.,
                             0.);

  /* Connect everything together */
  sys.addBus(&bus1);
  sys.addBus(&bus2);
  sys.addComponent(&branch);
  sys.addComponent(&fault);
  sys.addComponent(&gen);
  sys.allocate();

  double dt = 1.0 / 4.0 / 60.0;

  struct OutputData
  {
    double ti, Vr, Vi, dw;
  };

  std::vector<OutputData> output;

  auto buffer_write_cb = [&](double t)
  {
    std::vector<double>& yval = sys.y();

    output.push_back({.ti = t,
                      .Vr = yval[0],
                      .Vi = yval[1],
                      .dw = yval[3]});
  };

  /* Set up simulation */
  Ida<double, size_t> ida(&sys);
  ida.configureSimulation();

  /* Run simulation */
  double start = static_cast<double>(clock());
  // ida.printOutputF(0, 0, buffer);
  ida.initializeSimulation(0.0, false);
  ida.runSimulation(1.0, std::round((1.0 - 0.0) / dt), buffer_write_cb);
  fault.setStatus(1);
  ida.initializeSimulation(1.0, false);
  ida.runSimulation(1.1, std::round((1.1 - 1.0) / dt), buffer_write_cb);
  fault.setStatus(0);
  ida.initializeSimulation(1.1, false);
  ida.runSimulation(10.0, std::round((10.0 - 1.1) / dt), buffer_write_cb);
  double stop = static_cast<double>(clock());

  double error_V = 0.0; // error in |V|

  // Read through the simulation data storred in the buffer
  for (size_t i = 0; i < output.size(); i++)
  {
    OutputData           data    = output[i];
    std::vector<double>& ref_sol = reference_solution[i + 1];

    double err =
        std::abs(std::sqrt(data.Vr * data.Vr + data.Vi * data.Vi) - ref_sol[2])
        / (1.0 + std::abs(ref_sol[2]));
    if (err > error_V)
      error_V = err;

    std::cout << "GridKit: t = " << data.ti
              << ", |V| = " << std::sqrt(data.Vr * data.Vr + data.Vi * data.Vi)
              << ", w = " << (1.0 + data.dw) << "\n";
    std::cout << "Ref    : t = " << ref_sol[0]
              << ", |V| = " << ref_sol[2]
              << ", w = " << ref_sol[1]
              << "\n";
    std::cout << "Error in |V| = "
              << err
              << "\n";
    std::cout << "\n";
  }

  int status = 0;
  std::cout << "Max error in |V| = " << error_V << "\n";
  if (error_V > 2e-4)
  {
    std::cout << "Test failed with error too large!\n";
    status = 1;
  }

  std::cout << "\n\nComplete in " << (stop - start) / CLOCKS_PER_SEC << " seconds\n";

  return status;
}
