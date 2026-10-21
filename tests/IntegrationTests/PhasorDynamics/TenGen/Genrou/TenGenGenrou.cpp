/**
 * @file TenGenGenrou.cpp
 * @author Adam Birchfield (abirchfield@tamu.edu)
 * @author Slaven Peles (peless@ornl.gov)
 * @brief Example running a 10-bus system
 *
 * Simulates a 11-bus system with 10 Genrou 6th order generator models
 * split in two groups of five generators. The two groups are connected
 * by a high-impedance branch, which makes connection between them weak.
 *
 */
#include <cmath>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <vector>

#include <GridKit/Model/PhasorDynamics/ComponentLibrary.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

using scalar_type = double;
using real_type   = double;
using index_type  = size_t;

/// Machine data shared by the example generators
GridKit::PhasorDynamics::GenrouData<real_type, index_type> genData(real_type p0, real_type q0)
{
  using Parameter = GridKit::PhasorDynamics::GenrouParameters;

  GridKit::PhasorDynamics::GenrouData<real_type, index_type> data;
  data.parameters[Parameter::p0] = p0;
  data.parameters[Parameter::q0] = q0;

  return data;
}

/// Branch data with the given line and transformer parameters
GridKit::PhasorDynamics::BranchData<real_type, index_type> branchData(const real_type R,
                                                                      const real_type X,
                                                                      const real_type G,
                                                                      const real_type B,
                                                                      const real_type tap   = 1.0,
                                                                      const real_type phase = 0.0)
{
  using Parameter = GridKit::PhasorDynamics::BranchParameters;

  GridKit::PhasorDynamics::BranchData<real_type, index_type> data;
  data.parameters[Parameter::R]     = R;
  data.parameters[Parameter::X]     = X;
  data.parameters[Parameter::G]     = G;
  data.parameters[Parameter::B]     = B;
  data.parameters[Parameter::tap]   = tap;
  data.parameters[Parameter::phase] = phase;

  return data;
}

int main()
{
  using namespace GridKit::PhasorDynamics;
  using namespace AnalysisManager::Sundials;
  using GridKit::Testing::isEqual;
  using Log = GridKit::Utilities::Logger;

  // Print summaries, such as the run time, without lowering a higher verbosity
  Log::raiseVerbosity(Log::SUMMARY);

  GridKit::Testing::TestStatus success = true;

  /* Create model parts */
  BusInfinite<scalar_type, index_type> bus1(1, 0);

  Bus<scalar_type, index_type> bus2(0.999733719609643, 0.0230757421433613);
  Bus<scalar_type, index_type> bus3(0.999047460127767, 0.0436368240395443);
  Bus<scalar_type, index_type> bus4(0.998097277498088, 0.061658775943942);
  Bus<scalar_type, index_type> bus5(0.997021447662007, 0.0771246581966912);
  Bus<scalar_type, index_type> bus6(0.946436885707683, 0.322888837484268);
  Bus<scalar_type, index_type> bus7(0.943037519659334, 0.332686393642377);
  Bus<scalar_type, index_type> bus8(0.940418229359708, 0.340019961013984);
  Bus<scalar_type, index_type> bus9(0.938638861502395, 0.344901620288291);
  Bus<scalar_type, index_type> bus10(0.937739191669114, 0.347340277548916);

  Branch<scalar_type, index_type> branch12(branchData(0.001, 0.005, 0, 0));
  Branch<scalar_type, index_type> branch23(branchData(0.001, 0.005, 0, 0));
  Branch<scalar_type, index_type> branch34(branchData(0.001, 0.005, 0, 0));
  Branch<scalar_type, index_type> branch45(branchData(0.001, 0.005, 0, 0));
  Branch<scalar_type, index_type> branch56(branchData(0.001, 0.1, 0, 0));
  Branch<scalar_type, index_type> branch67(branchData(0.001, 0.005, 0, 0));
  Branch<scalar_type, index_type> branch78(branchData(0.001, 0.005, 0, 0));
  Branch<scalar_type, index_type> branch89(branchData(0.001, 0.005, 0, 0));
  Branch<scalar_type, index_type> branch910(branchData(0.001, 0.005, 0, 0));
  branch12.getPorts().bus.port<BranchBuses::bus1>().connect(&bus1);
  branch12.getPorts().bus.port<BranchBuses::bus2>().connect(&bus2);
  branch23.getPorts().bus.port<BranchBuses::bus1>().connect(&bus2);
  branch23.getPorts().bus.port<BranchBuses::bus2>().connect(&bus3);
  branch34.getPorts().bus.port<BranchBuses::bus1>().connect(&bus3);
  branch34.getPorts().bus.port<BranchBuses::bus2>().connect(&bus4);
  branch45.getPorts().bus.port<BranchBuses::bus1>().connect(&bus4);
  branch45.getPorts().bus.port<BranchBuses::bus2>().connect(&bus5);
  branch56.getPorts().bus.port<BranchBuses::bus1>().connect(&bus5);
  branch56.getPorts().bus.port<BranchBuses::bus2>().connect(&bus6);
  branch67.getPorts().bus.port<BranchBuses::bus1>().connect(&bus6);
  branch67.getPorts().bus.port<BranchBuses::bus2>().connect(&bus7);
  branch78.getPorts().bus.port<BranchBuses::bus1>().connect(&bus7);
  branch78.getPorts().bus.port<BranchBuses::bus2>().connect(&bus8);
  branch89.getPorts().bus.port<BranchBuses::bus1>().connect(&bus8);
  branch89.getPorts().bus.port<BranchBuses::bus2>().connect(&bus9);
  branch910.getPorts().bus.port<BranchBuses::bus1>().connect(&bus9);
  branch910.getPorts().bus.port<BranchBuses::bus2>().connect(&bus10);

  Genrou<scalar_type, index_type> gen2(genData(0.5, -0.00442101));
  Genrou<scalar_type, index_type> gen3(genData(0.5, -0.02510812));
  Genrou<scalar_type, index_type> gen4(genData(0.5, -0.04339553));
  Genrou<scalar_type, index_type> gen5(genData(0.5, -0.2334993));
  Genrou<scalar_type, index_type> gen6(genData(0.5, 0.69907194));
  Genrou<scalar_type, index_type> gen7(genData(0.5, -0.08318208));
  Genrou<scalar_type, index_type> gen8(genData(0.5, -0.09123614));
  Genrou<scalar_type, index_type> gen9(genData(0.5, -0.09662372));
  Genrou<scalar_type, index_type> gen10(genData(0.5, -0.09932297));
  gen2.getPorts().bus.port<GenrouBuses::bus>().connect(&bus2);
  gen3.getPorts().bus.port<GenrouBuses::bus>().connect(&bus3);
  gen4.getPorts().bus.port<GenrouBuses::bus>().connect(&bus4);
  gen5.getPorts().bus.port<GenrouBuses::bus>().connect(&bus5);
  gen6.getPorts().bus.port<GenrouBuses::bus>().connect(&bus6);
  gen7.getPorts().bus.port<GenrouBuses::bus>().connect(&bus7);
  gen8.getPorts().bus.port<GenrouBuses::bus>().connect(&bus8);
  gen9.getPorts().bus.port<GenrouBuses::bus>().connect(&bus9);
  gen10.getPorts().bus.port<GenrouBuses::bus>().connect(&bus10);

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
  sys.allocate();

  real_type dt = 1.0 / 4.0 / 60.0;

  // Uncomment code below to print output to a file:
  std::ofstream fileout;
  fileout.open("TenGenGenrou_Results.csv");
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
    const auto* yval = sys.y().getData();

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
      // We are outputting second equation of each generator.
      out << yval[18 + gen2.size() * i + 1] << ",";
    }
    out << "\n";
  };

  // Set up simulation
  Ida<scalar_type, index_type> ida(&sys);
  ida.configureSimulation();

  // Run simulation, output each `dt` interval
  real_type start = static_cast<real_type>(clock());
  ida.initializeSimulation(0.0, false);

  const auto*                  y = sys.y().getData();
  const std::vector<real_type> initial(y, y + sys.y().getSize());

  // Run for 1s and check that the initial operating point is preserved.
  ida.runSimulation(1.0, dt, output_cb);
  y = sys.y().getData();
  for (index_type i = 0; i < sys.y().getSize(); ++i)
  {
    success *= isEqual(y[i], initial[i], 1e-4);
  }

  // Introduce fault to ground and run for 0.1s
  bus10.setFault(true, 0.0, 1e-5);
  ida.initializeSimulation(1.0);
  ida.runSimulation(1.1, dt, output_cb);

  // Pm = 0.5 pu, H = 3 s; neglect electrical output and damping during the fault.
  using Variable               = GenrouInternalVariables;
  const auto      omega_index  = static_cast<index_type>(Variable::OMEGA);
  const real_type omega_ref    = std::sqrt(1.0 + 0.5 * 0.1 / 3.0) - 1.0;
  success                     *= std::hypot(bus10.Vr(), bus10.Vi()) < 1e-3;
  success                     *= isEqual(gen10.y().getData()[omega_index], omega_ref, 5e-5);

  // Clear fault and run until t = 10s.
  bus10.setFault(false, 0.0, 1e-5);
  ida.initializeSimulation(1.1);
  ida.runSimulation(10.0, dt, output_cb);
  real_type stop = static_cast<real_type>(clock());

  // Voltage should recover to within 2% of its initial magnitude.
  const real_type voltage  = std::hypot(bus10.Vr(), bus10.Vi());
  success                 *= std::abs(voltage - 1.0) < 0.02;

  fileout.close();

  Log::summary() << "Complete in " << (stop - start) / CLOCKS_PER_SEC << " seconds\n";

  auto error_set = GridKit::Testing::compareCSV(
      "TenGenGenrou_Results.csv",
      "reference/TenGenGenrou.ref.csv",
      GridKit::Testing::ErrorType::ABSOLUTE);
  error_set->display();
  success *= error_set->total_error.max_value < 1e-4;

  return success.report("TenGenGenrou");
}
