/**
 * @file TenGenClassical.cpp
 * @author Adam Birchfield (abirchfield@tamu.edu)
 * @brief Example running a 10-bus system
 *
 * Simulates a 11-bus system with 10 GenClassical 2nd order generator models
 * split in two groups of five generators. The two groups are connected
 * by a high-impedance branch, which makes connection between them weak.
 * It is the same as Example_TenGen_Genrou but with simpler generator models
 *
 */
#include <cstdio>
#include <ctime>
#include <fstream>
#include <vector>

// #include "TenGenClassical.hpp"
#include <GridKit/Model/PhasorDynamics/ComponentLibrary.hpp>
#include <GridKit/Model/PhasorDynamics/SystemModel.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>
#include <GridKit/Testing/Testing.hpp>

using scalar_type = double;
using real_type   = double;
using index_type  = size_t;

int main()
{
  using namespace GridKit::PhasorDynamics;
  using namespace AnalysisManager::Sundials;

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

  Branch<scalar_type, index_type> branch12(0.001, 0.005, 0, 0);
  Branch<scalar_type, index_type> branch23(0.001, 0.005, 0, 0);
  Branch<scalar_type, index_type> branch34(0.001, 0.005, 0, 0);
  Branch<scalar_type, index_type> branch45(0.001, 0.005, 0, 0);
  Branch<scalar_type, index_type> branch56(0.001, 0.1, 0, 0);
  Branch<scalar_type, index_type> branch67(0.001, 0.005, 0, 0);
  Branch<scalar_type, index_type> branch78(0.001, 0.005, 0, 0);
  Branch<scalar_type, index_type> branch89(0.001, 0.005, 0, 0);
  Branch<scalar_type, index_type> branch910(0.001, 0.005, 0, 0);

  GenClassical<scalar_type, index_type> gen2(0.5, -0.00442101, 3., 0.1, 0., 0.2);
  GenClassical<scalar_type, index_type> gen3(0.5, -0.02510812, 3., 0.1, 0., 0.2);
  GenClassical<scalar_type, index_type> gen4(0.5, -0.04339553, 3., 0.1, 0., 0.2);
  GenClassical<scalar_type, index_type> gen5(0.5, -0.2334993, 3., 0.1, 0., 0.2);
  GenClassical<scalar_type, index_type> gen6(0.5, 0.69907194, 3., 0.1, 0., 0.2);
  GenClassical<scalar_type, index_type> gen7(0.5, -0.08318208, 3., 0.1, 0., 0.2);
  GenClassical<scalar_type, index_type> gen8(0.5, -0.09123614, 3., 0.1, 0., 0.2);
  GenClassical<scalar_type, index_type> gen9(0.5, -0.09662372, 3., 0.1, 0., 0.2);
  GenClassical<scalar_type, index_type> gen10(0.5, -0.09932297, 3., 0.1, 0., 0.2);

  BusFault<scalar_type, index_type> fault(0, 1e-5, 0);

  /* Wire components to the bus voltage signal nodes */
  SignalNode<scalar_type, index_type> bus_vr[10];
  SignalNode<scalar_type, index_type> bus_vi[10];

  BusBase<scalar_type, index_type>* buses[10] = {&bus1, &bus2, &bus3, &bus4, &bus5, &bus6, &bus7, &bus8, &bus9, &bus10};
  for (int i = 0; i < 10; ++i)
  {
    buses[i]->getSignals().assignSignalNode<BusInternalVariables::VR>(&bus_vr[i]);
    buses[i]->getSignals().assignSignalNode<BusInternalVariables::VI>(&bus_vi[i]);
  }

  auto wire_branch = [&](Branch<scalar_type, index_type>& branch, int i, int j)
  {
    branch.getSignals().attachSignalNode<BranchExternalVariables::VR1>(&bus_vr[i]);
    branch.getSignals().attachSignalNode<BranchExternalVariables::VI1>(&bus_vi[i]);
    branch.getSignals().attachSignalNode<BranchExternalVariables::VR2>(&bus_vr[j]);
    branch.getSignals().attachSignalNode<BranchExternalVariables::VI2>(&bus_vi[j]);
  };
  wire_branch(branch12, 0, 1);
  wire_branch(branch23, 1, 2);
  wire_branch(branch34, 2, 3);
  wire_branch(branch45, 3, 4);
  wire_branch(branch56, 4, 5);
  wire_branch(branch67, 5, 6);
  wire_branch(branch78, 6, 7);
  wire_branch(branch89, 7, 8);
  wire_branch(branch910, 8, 9);

  auto wire_gen = [&](GenClassical<scalar_type, index_type>& gen, int i)
  {
    gen.getSignals().attachSignalNode<GenClassicalExternalVariables::VR>(&bus_vr[i]);
    gen.getSignals().attachSignalNode<GenClassicalExternalVariables::VI>(&bus_vi[i]);
  };
  wire_gen(gen2, 1);
  wire_gen(gen3, 2);
  wire_gen(gen4, 3);
  wire_gen(gen5, 4);
  wire_gen(gen6, 5);
  wire_gen(gen7, 6);
  wire_gen(gen8, 7);
  wire_gen(gen9, 8);
  wire_gen(gen10, 9);

  fault.getSignals().attachSignalNode<BusFaultExternalVariables::VR>(&bus_vr[9]);
  fault.getSignals().attachSignalNode<BusFaultExternalVariables::VI>(&bus_vi[9]);

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
  fileout.open("TenGenClassical_results.csv");
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
  for (size_t i = 0; i < 9; ++i)
  {
    out << "delta" << i + 2 << ",";
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
      // Each generator has 5 variables.
      // We are outputting second variables of each generator.
      out << yval[18 + 5 * i + 1] << ",";
    }

    // Output generator angles
    for (size_t i = 0; i < 9; ++i)
    {
      // 18 is offset for variables of 9 buses.
      // Each generator has 5 variables.
      // We are outputting first variables of each generator.
      out << yval[18 + 5 * i] << ",";
    }
    out << "\n";
  };

  // Set up simulation
  Ida<scalar_type, index_type> ida(&sys);
  ida.configureSimulation();

  // Run simulation, output each `dt` interval
  real_type start = static_cast<real_type>(clock());
  ida.initializeSimulation(0.0, false);

  // Run for 1s
  ida.runSimulation(1.0, dt, output_cb);

  // Introduce fault to ground and run for 0.1s
  fault.setStatus(1);
  ida.initializeSimulation(1.0);
  ida.runSimulation(1.1, dt, output_cb);

  // Clear fault and run until t = 10s.
  fault.setStatus(0);
  ida.initializeSimulation(1.1);
  ida.runSimulation(10.0, dt, output_cb);
  real_type stop = static_cast<real_type>(clock());

  fileout.close();

  std::cout << "\n\nComplete in " << (stop - start) / CLOCKS_PER_SEC << " seconds\n";

  return 0;
}
