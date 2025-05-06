/**
 * @file example2.cpp
 * @author Adam Birchfield (abirchfield@tamu.edu)
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Example running a 3-bus system
 *
 * Simulates a 3-bus system with two Genrou 6th order generator models and
 * compares results with data generated for the same system by Poweworld.
 *
 */
#include "example2.hpp"

#include <cstdio>
#include <ctime>
#include <fstream>
#include <vector>

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

struct OutputData
{
  real_type   t;
  scalar_type gen2speed;
  scalar_type gen3speed;
  scalar_type v2mag;
  scalar_type v3mag;

  OutputData operator-(const OutputData& other) const
  {
    assert(t == other.t);
    return {
        .t         = t,
        .gen2speed = gen2speed - other.gen2speed,
        .gen3speed = gen3speed - other.gen3speed,
        .v2mag     = v2mag - other.v2mag,
        .v3mag     = v3mag - other.v3mag,
    };
  }

  double norm() const
  {
    return std::max({
        std::abs(gen2speed),
        std::abs(gen3speed),
        std::abs(v2mag),
        std::abs(v3mag),
    });
  }
};

std::ostream& operator<<(std::ostream& out, const OutputData& data)
{
  out << data.t << ","
      << data.gen2speed << ","
      << data.gen3speed << ","
      << data.v2mag << ","
      << data.v3mag;
  return out;
}

int main()
{
  using namespace GridKit::PhasorDynamics;
  using namespace AnalysisManager::Sundials;

  auto error_allowed = static_cast<real_type>(1e-4);

  std::cout << "Example 2 version 1\n";

  /* Create model parts */
  SystemModel<scalar_type, index_type> sys;
  BusInfinite<scalar_type, index_type> bus1(1.06, 0.0);
  Bus<scalar_type, index_type>         bus2(1.0599558398065716, -0.009675621941024773);
  Bus<scalar_type, index_type>         bus3(0.9610827543495831, -0.13122476630506485);
  Branch<scalar_type, index_type>      branch12(&bus1, &bus2, 0.05, 0.21, 0, 0.1);
  Branch<scalar_type, index_type>      branch13(&bus1, &bus3, 0.06, 0.15, 0, 0.12);
  Branch<scalar_type, index_type>      branch23(&bus2, &bus3, 0.08, 0.27, 0, 0.45);
  Genrou<scalar_type, index_type>      gen2(&bus2, 1, 0.5, -0.07588, 2.7, 0., 0., 7., .04, .05, .75, 1.9, 0.17, 0.15, 0.4, 0.35, 0.15, 0.14999, 0., 0.);
  Genrou<scalar_type, index_type>      gen3(&bus3, 1, 0.25, 0.26587, 1.6, 0., 0., 7.5, .04, .05, .75, 2.3, 0.2, 0.18, 0.5, 0.5, 0.18, 0.15, 0., 0.);
  Load<scalar_type, index_type>        load3(&bus3, 0.4447197839297772, 0.20330047265361242);
  BusFault<scalar_type, index_type>    fault(&bus3, 0, 1e-5, 0);

  /* Connect everything together */
  sys.addBus(&bus1);
  sys.addBus(&bus2);
  sys.addBus(&bus3);
  sys.addComponent(&branch12);
  sys.addComponent(&branch13);
  sys.addComponent(&branch23);
  sys.addComponent(&fault);
  sys.addComponent(&load3);
  sys.addComponent(&gen2);
  sys.addComponent(&gen3);
  sys.allocate();

  real_type dt = 1.0 / 4.0 / 60.0;

  std::vector<OutputData> output;

  auto output_cb = [&](real_type t)
  {
    output.push_back({.t         = t,
                      .gen2speed = 1 + gen2.y()[1],
                      .gen3speed = 1 + gen3.y()[1],
                      .v2mag     = sqrt(bus2.Vr() * bus2.Vr() + bus2.Vi() * bus2.Vi()),
                      .v3mag     = sqrt(bus3.Vr() * bus3.Vr() + bus3.Vi() * bus3.Vi())});
  };

  /* Set up simulation */
  Ida<scalar_type, index_type>
      ida(&sys);
  ida.configureSimulation();

  /* Run simulation */
  scalar_type start = static_cast<scalar_type>(clock());
  ida.initializeSimulation(0.0, false);
  ida.runSimulation(1.0, std::round((1.0 - 0.0) / dt), output_cb);
  fault.setStatus(1);
  ida.initializeSimulation(1.0, false);
  ida.runSimulation(1.1, std::round((1.1 - 1.0) / dt), output_cb);
  fault.setStatus(0);
  ida.initializeSimulation(1.1, false);
  ida.runSimulation(10.0, std::round((10.0 - 1.1) / dt), output_cb);
  double stop = static_cast<double>(clock());

  /* Check worst-case error */
  real_type worst_error      = 0;
  real_type worst_error_time = 0;

  std::ostream  nullout(nullptr);
  std::ostream& out = nullout;

  // // Uncomment code below to print output to a file:
  // std::ofstream fileout;
  // fileout.open("example2_results.csv");
  // std::ostream& out = fileout;

  out << "Time,gen2speed,gen3speed,v2mag,v3mag\n";
  out << 0. << "," << 1. << "," << 1. << "," << 1. << "," << 1. << "\n";

  for (index_type i = 0; i < output.size(); ++i)
  {
    OutputData ref = {
        .t         = reference_solution[i + 1][0],
        .gen2speed = reference_solution[i + 1][1],
        .gen3speed = reference_solution[i + 1][2],
        .v2mag     = reference_solution[i + 1][4],
        .v3mag     = reference_solution[i + 1][5]};
    OutputData out_data = output[i];

    out << out_data << '\n';

    real_type err = (out_data - ref).norm();
    if (err > worst_error)
    {
      worst_error      = err;
      worst_error_time = out_data.t;
    }
  }
  // fileout.close();

  std::cout << "Worst error " << worst_error
            << " at time t = " << worst_error_time << "\n";
  std::cout << "\n\nComplete in " << (stop - start) / CLOCKS_PER_SEC << " seconds\n";

  return worst_error < error_allowed ? 0 : 1;
}
