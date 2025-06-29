/**
 * @file ThreeBusBasic.cpp
 * @author Adam Birchfield (abirchfield@tamu.edu)
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Example running a 3-bus system
 *
 * Simulates a 3-bus system with two Genrou 6th order generator models and
 * compares results with data generated for the same system by Poweworld.
 *
 */
#include <cstdio>
#include <ctime>
#include <fstream>
#include <vector>

#include "ThreeBusBasic.hpp"
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
#include <Model/PhasorDynamics/SystemModelData.hpp>
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
    assert(GridKit::Testing::isEqual(t, other.t, reference_tol));
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

  std::cout << "Example: ThreeBusBasic\n";

  //
  // Create model data
  //

  SystemModelData<scalar_type, index_type> data;

  // Set bus data
  data.bus.resize(3);

  // Bus 0
  data.bus[0].bus_id   = 0;
  data.bus[0].bus_type = BusData<scalar_type, index_type>::BusType::SLACK;
  data.bus[0].Vr0      = 1.06;
  data.bus[0].Vi0      = 0.0;

  // Bus 1
  data.bus[1].bus_id   = 1;
  data.bus[1].bus_type = BusData<scalar_type, index_type>::BusType::DEFAULT;
  data.bus[1].Vr0      = 1.0599558398065716;
  data.bus[1].Vi0      = -0.009675621941024773;

  // Bus 2
  data.bus[2].bus_id   = 2;
  data.bus[2].bus_type = BusData<scalar_type, index_type>::BusType::DEFAULT;
  data.bus[2].Vr0      = 0.9610827543495831;
  data.bus[2].Vi0      = -0.13122476630506485;

  // Set branch data
  data.branch.resize(3);

  // Branch 0-1
  data.branch[0].bus1_id = data.bus[0].bus_id;
  data.branch[0].bus2_id = data.bus[1].bus_id;
  data.branch[0].R       = 0.05;
  data.branch[0].X       = 0.21;
  data.branch[0].G       = 0;
  data.branch[0].B       = 0.1;

  // Branch 0-2
  data.branch[1].bus1_id = data.bus[0].bus_id;
  data.branch[1].bus2_id = data.bus[2].bus_id;
  data.branch[1].R       = 0.06;
  data.branch[1].X       = 0.15;
  data.branch[1].G       = 0;
  data.branch[1].B       = 0.12;

  // Branch 1-2
  data.branch[2].bus1_id = data.bus[1].bus_id;
  data.branch[2].bus2_id = data.bus[2].bus_id;
  data.branch[2].R       = 0.08;
  data.branch[2].X       = 0.27;
  data.branch[2].G       = 0;
  data.branch[2].B       = 0.45;

  // Set generator data
  data.genrou.resize(2);

  // Generator on bus 1
  data.genrou[0].bus_id = 1;
  data.genrou[0].p0     = 0.5;
  data.genrou[0].q0     = -0.07588;
  data.genrou[0].H      = 2.7;
  data.genrou[0].D      = 0.;
  data.genrou[0].Ra     = 0.;
  data.genrou[0].Tdop   = 7.;
  data.genrou[0].Tdopp  = .04;
  data.genrou[0].Tqopp  = .05;
  data.genrou[0].Tqop   = .75;
  data.genrou[0].Xd     = 1.9;
  data.genrou[0].Xdp    = 0.17;
  data.genrou[0].Xdpp   = 0.15;
  data.genrou[0].Xq     = 0.4;
  data.genrou[0].Xqp    = 0.35;
  data.genrou[0].Xqpp   = 0.15;
  data.genrou[0].Xl     = 0.14999;
  data.genrou[0].S10    = 0.;
  data.genrou[0].S12    = 0.;

  // Generator on bus 2
  data.genrou[1].bus_id = 2;
  data.genrou[1].p0     = 0.25;
  data.genrou[1].q0     = 0.26587;
  data.genrou[1].H      = 1.6;
  data.genrou[1].D      = 0.;
  data.genrou[1].Ra     = 0.;
  data.genrou[1].Tdop   = 7.5;
  data.genrou[1].Tdopp  = .04;
  data.genrou[1].Tqopp  = .05;
  data.genrou[1].Tqop   = .75;
  data.genrou[1].Xd     = 2.3;
  data.genrou[1].Xdp    = 0.2;
  data.genrou[1].Xdpp   = 0.18;
  data.genrou[1].Xq     = 0.5;
  data.genrou[1].Xqp    = 0.5;
  data.genrou[1].Xqpp   = 0.18;
  data.genrou[1].Xl     = 0.15;
  data.genrou[1].S10    = 0.;
  data.genrou[1].S12    = 0.;

  // Set load data
  data.load.resize(1);

  // Load on bus 2
  data.load[0].bus_id = 2;
  data.load[0].R      = 0.4447197839297772;
  data.load[0].X      = 0.20330047265361242;

  // Set fault data
  data.bus_fault.resize(1);

  data.bus_fault[0].bus_id = 2;
  data.bus_fault[0].R      = 0.0;
  data.bus_fault[0].X      = 1e-5;
  data.bus_fault[0].status = false;

  //
  // Instantiate system
  //

  SystemModel<scalar_type, index_type> sys(data);
  sys.allocate();

  // Get access to fault 0
  auto* fault = sys.getBusFault(0);

  real_type dt = 1.0 / 4.0 / 60.0;

  std::vector<OutputData> output;

  auto output_cb = [&](real_type t)
  {
    std::vector<double>& y_val = sys.y();

    output.push_back(OutputData{t,
                                1.0 + y_val[5],
                                1.0 + y_val[25],
                                std::sqrt(y_val[0] * y_val[0] + y_val[1] * y_val[1]),
                                std::sqrt(y_val[2] * y_val[2] + y_val[3] * y_val[3])});
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
  fault->setStatus(true);
  ida.initializeSimulation(1.0, false);
  nout = static_cast<int>(std::round((1.1 - 1.0) / dt));
  ida.runSimulation(1.1, nout, output_cb);

  // Clear fault and run until t = 10s.
  fault->setStatus(false);
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
  // fileout.open("Example_ThreeBus_Basic_results.csv");
  // std::ostream& out = fileout;

  out << "Time,gen2speed,gen3speed,v2mag,v3mag\n";
  out << 0. << "," << 1. << "," << 1. << "," << 1. << "," << 1. << "\n";

  for (index_type i = 0; i < output.size(); ++i)
  {
    OutputData ref{reference_solution[i + 1][0],
                   reference_solution[i + 1][1],
                   reference_solution[i + 1][2],
                   reference_solution[i + 1][4],
                   reference_solution[i + 1][5]};
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
