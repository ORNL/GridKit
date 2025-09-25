/**
 * @file TwoBusIeeet1.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @author Adam Birchfield (abirchfield@tamu.edu)
 * @brief Example running a 2-bus system with exciter and governor
 *
 */
#include "TwoBusIeeet1.hpp"

#include <ctime>
#include <iostream>

#include <Model/PhasorDynamics/ComponentLibrary.hpp>
#include <Model/PhasorDynamics/SystemModel.hpp>
#include <Model/PhasorDynamics/SystemModelData.hpp>
#include <Solver/Dynamic/Ida.hpp>
#include <Utilities/Testing.hpp>

// Temp, remove
#include <Model/PhasorDynamics/Exciter/IEEET1/Ieeet1Data.hpp>

int main()
{
  using namespace GridKit::PhasorDynamics;
  using namespace AnalysisManager::Sundials;
  using namespace GridKit::PhasorDynamics::Governor;
  using namespace GridKit::PhasorDynamics::Exciter;

  using scalar_type = double;
  using real_type   = double;
  using index_type  = size_t;

  using BusType      = BusData<scalar_type, index_type>::BusType;
  using signal_type  = SignalNode<scalar_type, index_type>;
  using machine_type = Genrou<scalar_type, index_type>;
  using gov_type     = Governor::Tgov1<scalar_type, index_type>;
  using exc_type     = Exciter::Ieeet1<scalar_type, index_type>;

  std::cout << "Example: TwoBusTgov1 + IEEET1 Exciter\n";

  // Model Data
  SystemModelData<scalar_type, index_type> data;

  // Set bus data
  data.bus.resize(2);

  data.bus[0].bus_id   = 0;
  data.bus[0].bus_type = BusType::DEFAULT;
  data.bus[0].Vr0      = 0.9949877346411762;
  data.bus[0].Vi0      = 0.09999703952427966;

  data.bus[1].bus_id   = 1;
  data.bus[1].bus_type = BusType::SLACK;
  data.bus[1].Vr0      = 1.0;
  data.bus[1].Vi0      = 0.0;

  // Set signal nodes data
  data.signal.resize(3);

  data.signal[0].name      = "omega";
  data.signal[0].signal_id = 0;

  data.signal[1].name      = "Pm";
  data.signal[1].signal_id = 1;

  data.signal[2].name      = "Efd";
  data.signal[2].signal_id = 2;

  // Set branch data
  data.branch.resize(1);

  data.branch[0].ports[BranchPorts::bus1]        = data.bus[0].bus_id;
  data.branch[0].ports[BranchPorts::bus2]        = data.bus[1].bus_id;
  data.branch[0].parameters[BranchParameters::R] = 0.0;
  data.branch[0].parameters[BranchParameters::X] = 0.1;
  data.branch[0].parameters[BranchParameters::G] = 0.0;
  data.branch[0].parameters[BranchParameters::B] = 0.0;

  // Add faults
  data.bus_fault.resize(1);

  data.bus_fault[0].parameters[BusFaultParameters::R]      = 0.0;
  data.bus_fault[0].parameters[BusFaultParameters::X]      = 1e-3;
  data.bus_fault[0].parameters[BusFaultParameters::state0] = false;

  // ------------- MODEL GROUP ------------------

  // Set generator data
  data.genrou.resize(1);

  data.genrou[0].ports[GenrouPorts::speed]           = 0;
  data.genrou[0].ports[GenrouPorts::pmech]           = 1;
  data.genrou[0].ports[GenrouPorts::efd]             = 2;
  data.genrou[0].parameters[GenrouParameters::p0]    = 1.;
  data.genrou[0].parameters[GenrouParameters::q0]    = 0.05013;
  data.genrou[0].parameters[GenrouParameters::H]     = 3.;
  data.genrou[0].parameters[GenrouParameters::D]     = 0.;
  data.genrou[0].parameters[GenrouParameters::Ra]    = 0.;
  data.genrou[0].parameters[GenrouParameters::Tdop]  = 7.;
  data.genrou[0].parameters[GenrouParameters::Tdopp] = .04;
  data.genrou[0].parameters[GenrouParameters::Tqopp] = .05;
  data.genrou[0].parameters[GenrouParameters::Tqop]  = .75;
  data.genrou[0].parameters[GenrouParameters::Xd]    = 2.1;
  data.genrou[0].parameters[GenrouParameters::Xdp]   = 0.2;
  data.genrou[0].parameters[GenrouParameters::Xdpp]  = 0.18;
  data.genrou[0].parameters[GenrouParameters::Xq]    = 0.5;
  data.genrou[0].parameters[GenrouParameters::Xqp]   = 0.5;
  data.genrou[0].parameters[GenrouParameters::Xqpp]  = 0.18;
  data.genrou[0].parameters[GenrouParameters::Xl]    = 0.15;
  data.genrou[0].parameters[GenrouParameters::S10]   = 0.;
  data.genrou[0].parameters[GenrouParameters::S12]   = 0.;

  // Governor
  data.gov.resize(1);

  data.gov[0].ports[Tgov1Ports::speed]           = 0;
  data.gov[0].ports[Tgov1Ports::pmech]           = 1;
  data.gov[0].parameters[Tgov1Parameters::R]     = 0.05;
  data.gov[0].parameters[Tgov1Parameters::Pvmin] = 0.0;
  data.gov[0].parameters[Tgov1Parameters::Pvmax] = 1.0;
  data.gov[0].parameters[Tgov1Parameters::T1]    = 0.5;
  data.gov[0].parameters[Tgov1Parameters::T2]    = 2.5;
  data.gov[0].parameters[Tgov1Parameters::T3]    = 7.5;
  data.gov[0].parameters[Tgov1Parameters::Dt]    = 0.0;

  // Exciter
  data.exciter.resize(1);

  data.exciter[0].ports[Ieeet1Ports::speed]                      = 0;
  data.exciter[0].ports[Ieeet1Ports::efd]                        = 2;
  data.exciter[0].parameters[Exciter::Ieeet1Parameters::Tr]      = 0.001; // (BUG: Nonfunctional if Tr = 0)
  data.exciter[0].parameters[Exciter::Ieeet1Parameters::Ka]      = 50.;
  data.exciter[0].parameters[Exciter::Ieeet1Parameters::Ta]      = 0.04;
  data.exciter[0].parameters[Exciter::Ieeet1Parameters::Ke]      = -0.06;
  data.exciter[0].parameters[Exciter::Ieeet1Parameters::Te]      = 0.6;
  data.exciter[0].parameters[Exciter::Ieeet1Parameters::Kf]      = 0.09;
  data.exciter[0].parameters[Exciter::Ieeet1Parameters::Tf]      = 1.46;
  data.exciter[0].parameters[Exciter::Ieeet1Parameters::Vrmin]   = -1.;
  data.exciter[0].parameters[Exciter::Ieeet1Parameters::Vrmax]   = 1.;
  data.exciter[0].parameters[Exciter::Ieeet1Parameters::E1]      = 2.8;
  data.exciter[0].parameters[Exciter::Ieeet1Parameters::E2]      = 3.373;
  data.exciter[0].parameters[Exciter::Ieeet1Parameters::Se1]     = 0.04;
  data.exciter[0].parameters[Exciter::Ieeet1Parameters::Se2]     = 0.33;
  data.exciter[0].parameters[Exciter::Ieeet1Parameters::Ispdlim] = 0.;

  // -------------- END MODEL DATA ------------------

  // Manually add components
  // This is a workaround since signal connections are not implemented in parser

  // Create buses
  auto* bus0 = BusFactory<scalar_type, index_type>::create(data.bus[0]);
  auto* bus1 = BusFactory<scalar_type, index_type>::create(data.bus[1]);

  // Create signal nodes
  signal_type omega(data.signal[0]);
  signal_type pmech(data.signal[1]);
  signal_type efd(data.signal[2]);

  // Create branch
  Branch<scalar_type, index_type> branch(bus0, bus1, data.branch[0]);

  // Add bus fault to bus0
  BusFault<scalar_type, index_type> fault(bus0, data.bus_fault[0]);

  // Create generator
  machine_type gen(bus0,   // Bus
                   &omega, // Machine  Speed Signal
                   &pmech, // Governor Pmech Signal
                   &efd,   // Exciter  Efd   Signal
                   data.genrou[0]);

  // Create governor (w/ Pmech and Speed signals)
  gov_type gov(&pmech, &omega);

  // Create exciter (w/ Efd, speed, and bus signals)
  exc_type exc(&efd, &omega, bus0, data.exciter[0]);

  // Instantiate system model and add components to it
  SystemModel<scalar_type, index_type> sys;
  sys.addBus(bus0);
  sys.addBus(bus1);
  sys.addSignal(&omega);
  sys.addSignal(&pmech);
  sys.addSignal(&efd);
  sys.addComponent(&branch);
  sys.addComponent(&gen);
  sys.addComponent(&gov);
  sys.addComponent(&exc);
  sys.addFault(&fault);

  sys.allocate();

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
    real_type ti, Vr, Vi, dw, Pm, Efd;
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

    // Note Omega of gen is at state index 5! (Each added signal shifted by 1)
    // Bus              -> 2 States
    // Genrou           -> 19 States -> Start Idx 2
    // Gov              -> 3 States  -> Start Idx 21
    // Exc              -> 9 States  -> Start Idx 24
    output.push_back(OutputData{
        t,
        y_val[0],  // Bus Vr
        y_val[1],  // Bus Vi
        y_val[3],  // Gen Speed
        y_val[23], // Gov Pmech
        y_val[26], // Exc Efd
    });
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
  fault.setStatus(true);
  ida.initializeSimulation(1.0, false);
  nout = static_cast<int>(std::round((1.1 - 1.0) / dt));
  ida.runSimulation(1.1, nout, output_cb);

  // Clear the fault and run until t = 10s.
  fault.setStatus(false);
  ida.initializeSimulation(1.1, false);
  nout = static_cast<int>(std::round((10.0 - 1.1) / dt));
  ida.runSimulation(10.0, nout, output_cb);
  real_type stop = static_cast<real_type>(clock());

  real_type error_V   = 0.0; // error in |V|
  real_type error_w   = 0.0; // error in rotor speed
  real_type error_Efd = 0.0; // error in exciter Efd

  // Output Headers
  // std::cout << "Time\t|V|\tSpeed\tPm\tEfd\tVreg" << "\n";
  // ref_sol[1] -> Speed
  // ref_sol[2] -> |V|
  // ref_sol[3] -> Efd

  // Read through the simulation data stored in the buffer.
  // Since we captured by reference, output should be available
  // for us to read here, outside the callback.
  for (size_t i = 0; i < output.size(); i++)
  {
    OutputData              data    = output[i];
    std::vector<real_type>& ref_sol = reference_solution[i + 1];

    // Review Note: I believe the denominator should not have +1
    real_type err =
        std::abs(std::sqrt(data.Vr * data.Vr + data.Vi * data.Vi) - ref_sol[2])
        / (1.0 + std::abs(ref_sol[2]));
    if (err > error_V)
      error_V = err;

    // Review Note: I believe the denominator should not have +1
    err = std::abs(1.0 + data.dw - ref_sol[1]) / (1.0 + ref_sol[1]);
    if (err > error_w)
      error_w = err;

    // Exciter Error
    err = std::abs(data.Efd - ref_sol[3]) / (ref_sol[3]);
    if (err > error_Efd)
      error_Efd = err;

    // Optional output
    /*
    std::cout  << data.ti
               << " " << std::sqrt(data.Vr * data.Vr + data.Vi * data.Vi)
               << " " << (1.0 + data.dw)
               << " " << (data.Pm)
               << " " << (data.Efd)
               << " " << (data.VR)
               << " " << (data.ksat)
               << "\n";
    // std::cout << "Ref    : t = " << ref_sol[0]
    //           << ", |V| = " << ref_sol[2]
    //           << ", w = " << ref_sol[1]
    //           << "\n";
    // std::cout << "Error in |V| = "
    //           << err
    //           << "\n";
    // std::cout << "\n";
    */
  }

  // Errors allowed for agreement with Powerworld results
  real_type error_V_allowed   = 2e-4;
  real_type error_w_allowed   = 1e-4;
  real_type error_Efd_allowed = 1e-2;

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
  std::cout << "Max error in  Efd  = " << error_Efd << "\n";
  if (error_Efd > error_Efd_allowed)
  {
    std::cout << "Test failed with error too large!\n";
    status = 1;
  }

  std::cout << "\n\nComplete in " << (stop - start) / CLOCKS_PER_SEC << " seconds\n";

  // Free Bus pointers manually
  delete bus0;
  delete bus1;

  return status;
}
