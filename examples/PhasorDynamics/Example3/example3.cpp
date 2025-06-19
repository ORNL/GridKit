/**
 * @file example3.cpp
 * @author Adam Birchfield (abirchfield@tamu.edu)
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Example running a 10-bus system
 *
 * Simulates a 11-bus system with 10 Genrou 6th order generator models
 * split in two groups of five generators. The two groups are connected
 * by a high-impedance branch, which makes connection between them weak.
 *
 */
#include <cstdio>
#include <ctime>
#include <fstream>
#include <vector>

// #include "example3.hpp"
#include <Model/PhasorDynamics/Branch/Branch.hpp>
#include <Model/PhasorDynamics/Bus/Bus.hpp>
#include <Model/PhasorDynamics/Bus/BusInfinite.hpp>
#include <Model/PhasorDynamics/BusFault/BusFault.hpp>
#include <Model/PhasorDynamics/Load/Load.hpp>
#include <Model/PhasorDynamics/SynchronousMachine/GENROUwS/Genrou.hpp>
#include <Model/PhasorDynamics/SystemModel.hpp>
#include <Solver/Dynamic/Ida.hpp>
#include <Utilities/Testing.hpp>

using scalar_type = double;
using real_type   = double;
using index_type  = size_t;

int main()
{
  using namespace GridKit::PhasorDynamics;
  using namespace AnalysisManager::Sundials;

  /* Create model parts */
  SystemModel<double, size_t> sys;
  Bus<double, size_t>         bus1(0.9949877346411762, 0.09999703952427966);
  BusInfinite<double, size_t> bus2(1.0, 0.0);
  Branch<double, size_t>      branch(&bus1, &bus2, 0, 0.1, 0, 0);
  BusFault<double, size_t>    fault(&bus1, 0, 1e-3, 0);

  // Decleration
  Genrou<double, size_t>* gen;
  TurbineGov<double, size_t>* gov;

  // Instatiation
  gen = new Genrou<double, size_t>(&bus1,
                             1,
                             gov
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

  // Governor of Generator
  TurbineGov<double, size_t> turb(&gen);

  /* Connect everything together */
  SystemModel<scalar_type, index_type> sys;

  sys.addBus(&bus1);
  sys.addBus(&bus2);
  sys.addBus(&bus3);
  sys.addBus(&bus4);
  sys.addBus(&bus5);
  sys.addBus(&bus6);
  sys.addBus(&bus7);
  sys.addBus(&bus8);
  sys.addBus(&bus9);
  sys.addBus(&bus10);
  sys.addComponent(&branch12);
  sys.addComponent(&branch23);
  sys.addComponent(&branch34);
  sys.addComponent(&branch45);
  sys.addComponent(&branch56);
  sys.addComponent(&branch67);
  sys.addComponent(&branch78);
  sys.addComponent(&branch89);
  sys.addComponent(&branch910);
  sys.addComponent(&gen2);
  sys.addComponent(&gen3);
  sys.addComponent(&gen4);
  sys.addComponent(&gen5);
  sys.addComponent(&gen6);
  sys.addComponent(&gen7);
  sys.addComponent(&gen8);
  sys.addComponent(&gen9);
  sys.addComponent(&gen10);
  sys.addComponent(&fault);
  sys.addComponent(&gen);
  sys.allocate();

  real_type dt = 1.0 / 4.0 / 60.0;

  // Uncomment code below to print output to a file:
  std::ofstream fileout;
  fileout.open("example3_results.csv");
  std::ostream& out = fileout;

  // Create header for the CSV output file
  out << "t,";
  for (size_t i = 0; i < 9; ++i)
  {
    out << "v" << i + 2 << ",";
  }
  for (size_t i = 0; i < 9; ++i)
  {
    out << "omega" << i + 2 << ",";
  }
  out << "\n";

  auto output_cb = [&](real_type t)
  {
    std::vector<double>& yval = sys.y();

    // Output time
    out << t << ",";

    // Output voltage magnitudes on buses
    for (size_t i = 0; i < 9; ++i)
    {
      out << std::sqrt(yval[2 * i] * yval[2 * i]
                       + yval[2 * i + 1] * yval[2 * i + 1])
          << ",";
    }

    // Output generator frequencies
    for (size_t i = 0; i < 9; ++i)
    {
      // 18 is offset for variables of 9 buses.
      // Each generator has 21 equations.
      // We are outputting second equation of each generator.
      out << yval[18 + 21 * i + 1] << ",";
    }
    out << "\n";
  };

  // Set up simulation
  Ida<scalar_type, index_type> ida(&sys);
  ida.configureSimulation();

  // Run simulation, output each `dt` interval
  scalar_type start = static_cast<scalar_type>(clock());
  ida.initializeSimulation(0.0, false);

  // Run for 1s
  int nout = static_cast<int>(std::round((1.0 - 0.0) / dt));
  ida.runSimulation(1.0, nout, output_cb);

  // Introduce fault to ground and run for 0.1s
  fault.setStatus(1);
  ida.initializeSimulation(1.0, false);
  nout = static_cast<int>(std::round((1.1 - 1.0) / dt));
  ida.runSimulation(1.1, nout, output_cb);

  // Clear fault and run until t = 10s.
  fault.setStatus(0);
  ida.initializeSimulation(1.1, false);
  nout = static_cast<int>(std::round((10.0 - 1.1) / dt));
  ida.runSimulation(10.0, nout, output_cb);
  double stop = static_cast<double>(clock());

  fileout.close();

  return 0;
}
