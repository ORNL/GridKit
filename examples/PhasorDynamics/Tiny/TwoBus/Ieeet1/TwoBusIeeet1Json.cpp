/**
 * @file TwoBusIeeet1.cpp
 * @author Luke Lowery (lukel@tamu.edu)
 * @author Adam Birchfield (abirchfield@tamu.edu)
 * @brief Example running a 2-bus system with exciter and governor
 *
 */
#include <ctime>
#include <iostream>

#include "TwoBusIeeet1.hpp"
#include <Model/PhasorDynamics/ComponentLibrary.hpp>
#include <Model/PhasorDynamics/SystemModel.hpp>
#include <Model/PhasorDynamics/SystemModelData.hpp>
#include <Model/PhasorDynamics/SystemModelDataJSONParser.hpp>
#include <Solver/Dynamic/Ida.hpp>
#include <Utilities/Testing.hpp>
#include <nlohmann/json.hpp>

int main()
{
  using namespace GridKit::PhasorDynamics;
  using namespace AnalysisManager::Sundials;
  using namespace GridKit::PhasorDynamics::Governor;
  using namespace GridKit::PhasorDynamics::Exciter;

  using scalar_type = double;
  using real_type   = double;
  using index_type  = size_t;

  std::cout << "Example: TwoBusTgov1 + IEEET1 Exciter\n";

  //
  // Input file
  //

  const char input_file[] =
      R"({
            "header": {
                "format_version": 0,
                "format_revision": 1,
                "case_name": "2-bus test case with governor and exciter",
                "case_description": "A two-bus test case for demonstrating the dynamics format",
                "case_comments": "This case is set up to monitor the voltage at both buses and the machine angle and speed",
                "freq_base": 60.0,
                "va_base": 100e6
            },
            "buses": [
                { 
                  "number": 0,
                  "class": "bus",
                  "name": "Bus 1",
                  "init": {
                            "Vr":0.9949877346411762,
                            "Vi":0.09999703952427966
                          },
                  "v_base": 115e3,
                  "mon": ["Vr", "Vi"]
                },
                { 
                  "number": 1,
                  "class": "infinite_bus",
                  "name": "Bus 2",
                  "init": {
                            "Vr":1.0,
                            "Vi":0.0
                          },
                  "v_base": 115e3
                }
            ],
            "signals": [
                { "signal_id": 0, "name": "Machine Speed Deviation"},
                { "signal_id": 1, "name": "Mechanical Power"},
                { "signal_id": 2, "name": "Excitation Field"}
            ],
            "devices": [
                { 
                  "class": "Branch",
                  "ports": {"bus1":0, "bus2":1},
                  "id": "BR1",
                  "params": {"R":0.0, "X":0.1, "G":0.0, "B":0.0}
                },
                {
                  "class": "Genrou",
                  "ports": {"bus":0, "speed": 0, "pmech":1, "efd": 2},
                  "id": "DV1",
                  "params": {"p0":1.0, "q0":0.05013, "H":3.0, "D":0.0, "Ra":0.0, "Tdop":7.0, "Tdopp":0.04, "Tqopp":0.05, "Tqop":0.75, "Xd":2.1, "Xdp":0.2, "Xdpp":0.18, "Xq":0.5, "Xqp": 0.5, "Xqpp":0.18, "Xl":0.15, "S10":0.0, "S12":0.0},
                  "mon": ["delta", "omega"]
                },
                {
                  "class": "Tgov1",
                  "ports": {"bus":0, "speed": 0, "pmech":1},
                  "id": "DV2",
                  "params": {"R":0.05, "T1":0.5, "T2":2.5, "T3":7.5, "Pvmin":0.0, "Pvmax":1.0, "Dt":0.0}
                },
                {
                  "class": "Ieeet1",
                  "ports": {"bus":0, "speed": 0, "efd":2},
                  "id": "DV3",
                  "params": {"Tr":0.001, "Ka":50.0, "Ta":0.04, "Ke":-0.06, "Te":0.6, "Kf":0.09, "Tf":1.46, "Vrmin":-1.0, "Vrmax":1.0, "E1":2.8, "E2":3.373, "Se1":0.04, "Se2":0.33, "Ispdlim":0.0}
                },
                { 
                  "class": "bus_fault",
                  "ports": {"bus":0},
                  "id": "0",
                  "params": {"state0": false, "R":0.0, "X":1e-3}
                }
            ]
        })";

  //
  // Create model data
  //

  SystemModelData<scalar_type, index_type> data = json::parse(input_file);

  //
  // Instantiate and configure the system model
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
  // fault.setStatus(true);
  fault->setStatus(true);
  ida.initializeSimulation(1.0, false);
  nout = static_cast<int>(std::round((1.1 - 1.0) / dt));
  ida.runSimulation(1.1, nout, output_cb);

  // Clear the fault and run until t = 10s.
  // fault.setStatus(false);
  fault->setStatus(false);
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

  return status;
}
