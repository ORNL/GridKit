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
#include <Model/PhasorDynamics/Branch/BranchData.hpp>
#include <Model/PhasorDynamics/Bus/Bus.hpp>
#include <Model/PhasorDynamics/Bus/BusData.hpp>
#include <Model/PhasorDynamics/Bus/BusInfinite.hpp>
#include <Model/PhasorDynamics/BusFault/BusFault.hpp>
#include <Model/PhasorDynamics/BusFault/BusFaultData.hpp>
#include <Model/PhasorDynamics/Load/Load.hpp>
#include <Model/PhasorDynamics/Load/LoadData.hpp>
#include <Model/PhasorDynamics/SynchronousMachine/GENROUwS/Genrou.hpp>
#include <Model/PhasorDynamics/SynchronousMachine/GENROUwS/GenrouData.hpp>
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

  OutputData& operator-=(const OutputData& other)
  {
    assert(GridKit::Testing::isEqual(t, other.t, Example2::reference_tol));
    gen2speed -= other.gen2speed;
    gen3speed -= other.gen3speed;
    v2mag     -= other.v2mag;
    v3mag     -= other.v3mag;
    return *this;
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

const OutputData operator-(const OutputData& lhs, const OutputData& rhs)
{
  return OutputData(lhs) -= rhs;
}

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

  //
  // Create (load) model data
  //

  // Bus 1
  BusData<real_type, index_type> bus_data_1;
  bus_data_1.Vr0 = 1.06;
  bus_data_1.Vi0 = 0.0;

  // Bus 2
  BusData<real_type, index_type> bus_data_2;
  bus_data_2.Vr0 = 1.0599558398065716;
  bus_data_2.Vi0 = -0.009675621941024773;

  // Bus 3
  BusData<real_type, index_type> bus_data_3;
  bus_data_3.Vr0 = 0.9610827543495831;
  bus_data_3.Vi0 = -0.13122476630506485;

  // Branch 1-2
  BranchData<real_type, index_type> branch_data_1_2;
  branch_data_1_2.R = 0.05;
  branch_data_1_2.X = 0.21;
  branch_data_1_2.G = 0;
  branch_data_1_2.B = 0.1;

  // Branch 1-3
  BranchData<real_type, index_type> branch_data_1_3;
  branch_data_1_3.R = 0.06;
  branch_data_1_3.X = 0.15;
  branch_data_1_3.G = 0;
  branch_data_1_3.B = 0.12;

  // Branch 2-3
  BranchData<real_type, index_type> branch_data_2_3;
  branch_data_2_3.R = 0.08;
  branch_data_2_3.X = 0.27;
  branch_data_2_3.G = 0;
  branch_data_2_3.B = 0.45;

  // Generator on bus 2
  GenrouData<real_type, index_type> gen_data_2;
  gen_data_2.unit_id = 1;
  gen_data_2.p0      = 0.5;
  gen_data_2.q0      = -0.07588;
  gen_data_2.H       = 2.7;
  gen_data_2.D       = 0.;
  gen_data_2.Ra      = 0.;
  gen_data_2.Tdop    = 7.;
  gen_data_2.Tdopp   = .04;
  gen_data_2.Tqopp   = .05;
  gen_data_2.Tqop    = .75;
  gen_data_2.Xd      = 1.9;
  gen_data_2.Xdp     = 0.17;
  gen_data_2.Xdpp    = 0.15;
  gen_data_2.Xq      = 0.4;
  gen_data_2.Xqp     = 0.35;
  gen_data_2.Xqpp    = 0.15;
  gen_data_2.Xl      = 0.14999;
  gen_data_2.S10     = 0.;
  gen_data_2.S12     = 0.;

  // Generator on bus 3
  GenrouData<real_type, index_type> gen_data_3;
  gen_data_3.unit_id = 1;
  gen_data_3.p0      = 0.25;
  gen_data_3.q0      = 0.26587;
  gen_data_3.H       = 1.6;
  gen_data_3.D       = 0.;
  gen_data_3.Ra      = 0.;
  gen_data_3.Tdop    = 7.5;
  gen_data_3.Tdopp   = .04;
  gen_data_3.Tqopp   = .05;
  gen_data_3.Tqop    = .75;
  gen_data_3.Xd      = 2.3;
  gen_data_3.Xdp     = 0.2;
  gen_data_3.Xdpp    = 0.18;
  gen_data_3.Xq      = 0.5;
  gen_data_3.Xqp     = 0.5;
  gen_data_3.Xqpp    = 0.18;
  gen_data_3.Xl      = 0.15;
  gen_data_3.S10     = 0.;
  gen_data_3.S12     = 0.;

  // Load on bus 3
  LoadData<real_type, index_type> load_data_3;
  load_data_3.R      = 0.4447197839297772;
  load_data_3.X      = 0.20330047265361242;
  load_data_3.bus_id = 3;

  BusFaultData<real_type, index_type> bus_fault_data_3;
  bus_fault_data_3.R = 0.0;
  bus_fault_data_3.X = 1e-5;
  bus_fault_data_3.status = 0;

  //
  // Instantiate model components
  //

  BusInfinite<scalar_type, index_type> bus1(bus_data_1);
  Bus<scalar_type, index_type>         bus2(bus_data_2);
  Bus<scalar_type, index_type>         bus3(bus_data_3);

  Branch<scalar_type, index_type> branch12(&bus1, &bus2, branch_data_1_2);
  Branch<scalar_type, index_type> branch13(&bus1, &bus3, branch_data_1_3);
  Branch<scalar_type, index_type> branch23(&bus2, &bus3, branch_data_2_3);

  Genrou<scalar_type, index_type>   gen2(&bus2, gen_data_2);
  Genrou<scalar_type, index_type>   gen3(&bus3, gen_data_3);
  Load<scalar_type, index_type>     load3(&bus3, load_data_3);
  BusFault<scalar_type, index_type> fault(&bus3, bus_fault_data_3);

  /* Connect everything together */
  SystemModel<scalar_type, index_type> sys;
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
    std::vector<double>& yval = sys.y();

    output.push_back(OutputData{t,
                                1.0 + yval[5],
                                1.0 + yval[26],
                                std::sqrt(yval[0] * yval[0] + yval[1] * yval[1]),
                                std::sqrt(yval[2] * yval[2] + yval[3] * yval[3])});
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
    OutputData ref{Example2::reference_solution[i + 1][0],
                   Example2::reference_solution[i + 1][1],
                   Example2::reference_solution[i + 1][2],
                   Example2::reference_solution[i + 1][4],
                   Example2::reference_solution[i + 1][5]};
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

  std::cout << "Max error " << worst_error
            << " at time t = " << worst_error_time << "\n";
  std::cout << "\n\nComplete in " << (stop - start) / CLOCKS_PER_SEC << " seconds\n";

  return worst_error < error_allowed ? 0 : 1;
}
