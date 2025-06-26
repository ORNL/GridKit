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
#include <Model/PhasorDynamics/Bus/BusElectric/BusInfinite.hpp>
#include <Model/PhasorDynamics/Bus/BusElectric/BusNetwork.hpp>
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
  BusInfinite<scalar_type, index_type> bus1(1, 0);

  BusNetwork<scalar_type, index_type> bus2(0.999733719609643, 0.0230757421433613);
  BusNetwork<scalar_type, index_type> bus3(0.999047460127767, 0.0436368240395443);
  BusNetwork<scalar_type, index_type> bus4(0.998097277498088, 0.061658775943942);
  BusNetwork<scalar_type, index_type> bus5(0.997021447662007, 0.0771246581966912);
  BusNetwork<scalar_type, index_type> bus6(0.946436885707683, 0.322888837484268);
  BusNetwork<scalar_type, index_type> bus7(0.943037519659334, 0.332686393642377);
  BusNetwork<scalar_type, index_type> bus8(0.940418229359708, 0.340019961013984);
  BusNetwork<scalar_type, index_type> bus9(0.938638861502395, 0.344901620288291);
  BusNetwork<scalar_type, index_type> bus10(0.937739191669114, 0.347340277548916);

  Branch<scalar_type, index_type> branch12(&bus1, &bus2, 0.001, 0.005, 0, 0);
  Branch<scalar_type, index_type> branch23(&bus2, &bus3, 0.001, 0.005, 0, 0);
  Branch<scalar_type, index_type> branch34(&bus3, &bus4, 0.001, 0.005, 0, 0);
  Branch<scalar_type, index_type> branch45(&bus4, &bus5, 0.001, 0.005, 0, 0);
  Branch<scalar_type, index_type> branch56(&bus5, &bus6, 0.001, 0.1, 0, 0);
  Branch<scalar_type, index_type> branch67(&bus6, &bus7, 0.001, 0.005, 0, 0);
  Branch<scalar_type, index_type> branch78(&bus7, &bus8, 0.001, 0.005, 0, 0);
  Branch<scalar_type, index_type> branch89(&bus8, &bus9, 0.001, 0.005, 0, 0);
  Branch<scalar_type, index_type> branch910(&bus9, &bus10, 0.001, 0.005, 0, 0);

  Genrou<scalar_type, index_type> gen2(&bus2, 1, 0.5, -0.00442101, 3., 0., 0., 7., .04, .05, .75, 2.1, 0.2, 0.18, 0.5, 0.5, 0.18, 0.15, 0., 0.);
  Genrou<scalar_type, index_type> gen3(&bus3, 1, 0.5, -0.02510812, 3., 0., 0., 7., .04, .05, .75, 2.1, 0.2, 0.18, 0.5, 0.5, 0.18, 0.15, 0., 0.);
  Genrou<scalar_type, index_type> gen4(&bus4, 1, 0.5, -0.04339553, 3., 0., 0., 7., .04, .05, .75, 2.1, 0.2, 0.18, 0.5, 0.5, 0.18, 0.15, 0., 0.);
  Genrou<scalar_type, index_type> gen5(&bus5, 1, 0.5, -0.2334993, 3., 0., 0., 7., .04, .05, .75, 2.1, 0.2, 0.18, 0.5, 0.5, 0.18, 0.15, 0., 0.);
  Genrou<scalar_type, index_type> gen6(&bus6, 1, 0.5, 0.69907194, 3., 0., 0., 7., .04, .05, .75, 2.1, 0.2, 0.18, 0.5, 0.5, 0.18, 0.15, 0., 0.);
  Genrou<scalar_type, index_type> gen7(&bus7, 1, 0.5, -0.08318208, 3., 0., 0., 7., .04, .05, .75, 2.1, 0.2, 0.18, 0.5, 0.5, 0.18, 0.15, 0., 0.);
  Genrou<scalar_type, index_type> gen8(&bus8, 1, 0.5, -0.09123614, 3., 0., 0., 7., .04, .05, .75, 2.1, 0.2, 0.18, 0.5, 0.5, 0.18, 0.15, 0., 0.);
  Genrou<scalar_type, index_type> gen9(&bus9, 1, 0.5, -0.09662372, 3., 0., 0., 7., .04, .05, .75, 2.1, 0.2, 0.18, 0.5, 0.5, 0.18, 0.15, 0., 0.);
  Genrou<scalar_type, index_type> gen10(&bus10, 1, 0.5, -0.09932297, 3., 0., 0., 7., .04, .05, .75, 2.1, 0.2, 0.18, 0.5, 0.5, 0.18, 0.15, 0., 0.);

  BusFault<scalar_type, index_type> fault(&bus10, 0, 1e-5, 0);

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

  std::cout << "\n\nComplete in " << (stop - start) / CLOCKS_PER_SEC << " seconds\n";

  return 0;
}
