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
#include <Model/PhasorDynamics/Branch/BranchData.hpp>
#include <Model/PhasorDynamics/Bus/Bus.hpp>
#include <Model/PhasorDynamics/Bus/BusData.hpp>
#include <Model/PhasorDynamics/Bus/BusInfinite.hpp>
#include <Model/PhasorDynamics/BusFault/BusFault.hpp>
#include <Model/PhasorDynamics/BusFault/BusFaultData.hpp>
#include <Model/PhasorDynamics/SynchronousMachine/GENROUwS/Genrou.hpp>
#include <Model/PhasorDynamics/SynchronousMachine/GENROUwS/GenrouData.hpp>
#include <Model/PhasorDynamics/SystemModel.hpp>
#include <Model/PhasorDynamics/SystemModelData.hpp>
#include <Solver/Dynamic/Ida.hpp>
#include <Utilities/Testing.hpp>

int main()
{
  using namespace GridKit::PhasorDynamics;
  using namespace AnalysisManager::Sundials;

  using scalar_type = double;
  using real_type   = double;
  using index_type  = size_t;

  std::cout << "Example 1 version 2\n";

  //
  // Create model data
  //

  SystemModelData<scalar_type, index_type> data;

  // Set bus data
  data.bus.resize(2);

  data.bus[0].bus_id   = 0;
  data.bus[0].bus_type = BusData<scalar_type, index_type>::DEFAULT;
  data.bus[0].Vr0      = 0.9949877346411762;
  data.bus[0].Vi0      = 0.09999703952427966;

  data.bus[1].bus_id   = 1;
  data.bus[1].bus_type = BusData<scalar_type, index_type>::SLACK;
  data.bus[1].Vr0      = 1.0;
  data.bus[1].Vi0      = 0.0;

  // Set branch data
  data.branch.resize(1);

  data.branch[0].bus1_id = data.bus[0].bus_id;
  data.branch[0].bus2_id = data.bus[1].bus_id;
  data.branch[0].R       = 0.0;
  data.branch[0].X       = 0.1;
  data.branch[0].G       = 0.0;
  data.branch[0].B       = 0.0;

  // Set generator data
  data.genrou.resize(1);

  data.genrou[0].unit_id = 1;
  data.genrou[0].p0      = 1.;
  data.genrou[0].q0      = 0.05013;
  data.genrou[0].H       = 3.;
  data.genrou[0].D       = 0.;
  data.genrou[0].Ra      = 0.;
  data.genrou[0].Tdop    = 7.;
  data.genrou[0].Tdopp   = .04;
  data.genrou[0].Tqopp   = .05;
  data.genrou[0].Tqop    = .75;
  data.genrou[0].Xd      = 2.1;
  data.genrou[0].Xdp     = 0.2;
  data.genrou[0].Xdpp    = 0.18;
  data.genrou[0].Xq      = 0.5;
  data.genrou[0].Xqp     = 0.5;
  data.genrou[0].Xqpp    = 0.18;
  data.genrou[0].Xl      = 0.15;
  data.genrou[0].S10     = 0.;
  data.genrou[0].S12     = 0.;

  // Add faults
  data.bus_fault.resize(1);

  data.bus_fault[0].R      = 0.0;
  data.bus_fault[0].X      = 1e-3;
  data.bus_fault[0].status = false;

  //
  // Instantiate system model
  //

  SystemModel<scalar_type, index_type> sys(data);
  sys.allocate();

  // Get access to the fault
  auto* fault = sys.getBusFault(0);

  // Set time step to 1/4 of a 60Hz cycle
  real_type dt = 1.0 / 4.0 / 60.0;

  // A data structure to keep track of the data we want to
  // compare to the reference solution. Rather than keeping
  // the entire solution vector at every time step around,
  // we instead narrow down exactly what we want to keep.
  //
  // Since this struct is "simple" enough (no constructors or
  // assignment operators, and "simple" members), it is a POD
  // (plain ol' data), which have some benefits in C++.
  struct OutputData
  {
    // Output variables are time, real and imaginary voltage and
    // frequency deviation
    real_type ti, Vr, Vi, dw;
  };

  // A list of output for each time step.
  std::vector<OutputData> output;

  // A callback which will be called by the integrator after
  // each time step. It will be told the time of the current
  // state, and it is allowed to access the up-to-date state
  // of the components, which are captured by a closure
  // due to the [&] notation (every variable that is referenced
  // by the callback that is external to the callback itself -
  // here output, bus1, and gen - will be considered a
  // reference to that variable inside the callback). We select
  // the subset of the output we're interested in recording and
  // push it into output, which is updated outside the callback.
  auto output_cb = [&](real_type t)
  {
    std::vector<scalar_type>& y_val = sys.y();

    output.push_back(OutputData{t, y_val[0], y_val[1], y_val[3]});
  };

  // Set up simulation
  Ida<scalar_type, size_t> ida(&sys);
  ida.configureSimulation();

  // Run simulation - making sure to pass the callback to record output
  real_type start = static_cast<real_type>(clock());

  // Run for 1s
  ida.initializeSimulation(0.0, false);
  int nout = static_cast<int>(std::round((1.0 - 0.0) / dt));
  ida.runSimulation(1.0, nout, output_cb);

  // Introduce fault and run for the next 0.1s
  fault->setStatus(true);
  ida.initializeSimulation(1.0, false);
  nout = static_cast<int>(std::round((1.1 - 1.0) / dt));
  ida.runSimulation(1.1, nout, output_cb);

  // Clear the fault and run until t = 10s.
  fault->setStatus(false);
  ida.initializeSimulation(1.1, false);
  nout = static_cast<int>(std::round((10.0 - 1.1) / dt));
  ida.runSimulation(10.0, nout, output_cb);
  real_type stop = static_cast<real_type>(clock());

  real_type error_V = 0.0; // error in |V|
  real_type error_w = 0.0; // error in rotor speed

  // Read through the simulation data stored in the buffer.
  // Since we captured by reference, output should be available
  // for us to read here, outside the callback.
  for (size_t i = 0; i < output.size(); i++)
  {
    OutputData              data    = output[i];
    std::vector<real_type>& ref_sol = Example1::reference_solution[i + 1];

    real_type err =
        std::abs(std::sqrt(data.Vr * data.Vr + data.Vi * data.Vi) - ref_sol[2])
        / (1.0 + std::abs(ref_sol[2]));
    if (err > error_V)
      error_V = err;

    err = std::abs(1.0 + data.dw - ref_sol[1]) / (1.0 + ref_sol[1]);
    if (err > error_w)
      error_w = err;

    // // Optional output
    // std::cout << "GridKit: t = " << data.ti
    //           << ", |V| = " << std::sqrt(data.Vr * data.Vr + data.Vi * data.Vi)
    //           << ", w = " << (1.0 + data.dw) << "\n";
    // std::cout << "Ref    : t = " << ref_sol[0]
    //           << ", |V| = " << ref_sol[2]
    //           << ", w = " << ref_sol[1]
    //           << "\n";
    // std::cout << "Error in |V| = "
    //           << err
    //           << "\n";
    // std::cout << "\n";
  }

  // Errors allowed for agreement with Powerworld results
  real_type error_V_allowed = 2e-4;
  real_type error_w_allowed = 1e-4;

  // Tolerances based on Powerworld reference accuracy
  int status = 0;
  std::cout << "Max error in |V| = " << error_V << "\n";
  if (error_V > error_V_allowed)
  {
    std::cout << "Test failed with error too large!\n";
    status = 1;
  }
  std::cout << "Max error in  w  = " << error_w << "\n";
  if (error_w > error_w_allowed)
  {
    std::cout << "Test failed with error too large!\n";
    status = 1;
  }

  std::cout << "\n\nComplete in " << (stop - start) / CLOCKS_PER_SEC << " seconds\n";

  return status;
}
