#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include <GridKit/Model/PowerElectronics/PartitionInterface/BusPartitionInterface.hpp>
#include <GridKit/Model/PowerElectronics/SubsystemModel.hpp>
#include <GridKit/Solver/Dynamic/DynamicSolver.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>
#include <GridKit/Solver/Dynamic/MultiStageCosimulation.hpp>

#include <examples/PowerElectronics/ExamplesHelper/JacTestHelper.hpp>
#include <examples/PowerElectronics/ExamplesHelper/MicrogridNetwork.hpp>
#include <examples/PowerElectronics/ExamplesHelper/PartitionUtilities.hpp>
#include <examples/PowerElectronics/ExamplesHelper/SystemAssembler.hpp>

using Component = GridKit::CircuitComponent<double, int>;
using Node      = GridKit::PowerElectronics::NodeBase<double, int>;
using Subsystem = GridKit::SubsystemModel<double, int>;
using System    = GridKit::PowerElectronicsModel<double, int>;

template <typename ScalarT>
std::vector<ScalarT> linspace(ScalarT start, ScalarT end, std::size_t num_points);

struct GenDataIndex
{
  bool refframe;
  int  idx_y11;
  int  idx_y12;
  int  idx_y13;
  int  idx_y14;
  int  idx_y15;
};

/*
 * Verify partitioned residual and Jacobian evaluation against the
 *        monolithic microgrid model.
 *
 * The microgrid is manually divided into two subsystems. A partition interface
 * is introduced at the boundary between bus 1 and the line connecting buses 1
 * and 2. Each subsystem is evaluated independently, and its residuals are gathered
 * and compared with the monolithic reference.
 */
int main()
{

  constexpr int    N_size           = 1000;
  constexpr size_t max_iterations   = 2;
  constexpr size_t num_stages       = 2;
  constexpr size_t num_intervals_1  = 50;
  constexpr size_t num_intervals_2  = 100;
  constexpr size_t num_partitions   = 128;
  constexpr size_t MAX_STATE_OUTPUT = 1299;
  bool             use_jac          = true;

  // ---------------------------------------------------------------------------
  // build the grid network and assemble the system model
  // ---------------------------------------------------------------------------
  ScaleMicrogridNetwork<double, int> network(N_size);

  auto* system = new System(use_jac);

  assembleSystem(network, *system);
  system->allocate();

  std::vector<GenDataIndex> genData(5);

  // Identify the location for the variables required to compute delta, P and Q for each generator
  // for (int i = 0; i < 4; i++)
  // {
  //   auto* gen     = network.generators[i];
  //   auto& gendata = genData[i];

  //   gendata.refframe = (i == 0);
  //   gendata.idx_y11  = gen->getNodeConnection(11);
  //   gendata.idx_y12  = gen->getNodeConnection(12);
  //   gendata.idx_y13  = gen->getNodeConnection(13);
  //   gendata.idx_y14  = gen->getNodeConnection(14);
  //   gendata.idx_y15  = gendata.refframe ? -1 : gen->getNodeConnection(15);
  // }

  std::cout << system->y().getSize() << std::endl;

  auto* y  = system->y().getData();
  auto* yp = system->yp().getData();

  // Create initial points for states
  for (int i = 0; i < system->size(); i++)
  {
    y[i]  = 0.0;
    yp[i] = 0.0;
  }

  // Create initial derivatives specifics generated in MATLAB
  for (int i = 0; i < 2 * N_size; i++)
  {
    const auto& params = network.DGParam_list[i];

    const int offset = (i == 0) ? 0 : 13 * i - 1;

    yp[offset + 2] = params.Vn_;
    yp[offset + 4] = params.Kpv_ * params.Vn_;
    yp[offset + 6] = (params.Kpc_ * params.Kpv_ * params.Vn_) / params.Lf_;
  }

  // since the intial P_com = 0, the set the intial vector to the reference frame
  y[network.dg_signal.getNodeConnection(0).idx_] = network.DGParam_list[0].wb_;

  system->y().setDataUpdated();
  system->yp().setDataUpdated();

  // Add fault to the system
  for (size_t i = 0; i < network.loads.size(); i += 2)
  {
    network.loads[i]->setResistance(3);
    network.loads[i]->setInductance(0.1);
  }

  system->updateTime(0.0, 0.0);
  system->evaluateResidual();

  const auto* residual = system->getResidual().getData();

  double max_residual = 0.0;
  int    max_index    = -1;

  for (int i = 0; i < system->size(); ++i)
  {
    const double value = std::abs(residual[i]);

    if (value > max_residual)
    {
      max_residual = value;
      max_index    = i;
    }
  }

  std::cout << "Initial system consistency:\n"
            << "  max|F| = " << max_residual << '\n'
            << "  index  = " << max_index << '\n';

  // --------------------------------------------------------------------------------
  // Create 2 Partitions and  a Partition interfaces to split to network
  // --------------------------------------------------------------------------------
  std::vector<Subsystem*> subsystems;
  partitionNetwork(network, subsystems, num_partitions);

  // ---------------------------------------------------------------------------
  // cosim
  // ---------------------------------------------------------------------------

  std::vector<double> Tspan1 = linspace(0.0, 0.01, num_intervals_1 + 1);
  std::vector<double> Tspan2 = linspace(0.01, 1.0, num_intervals_2 + 1);
  std::vector<double> Tspan  = Tspan1;

  Tspan.insert(Tspan.end(), Tspan2.begin() + 1, Tspan2.end());

  auto* cosim = new AnalysisManager::MultiStageCosimulation<double, int>(max_iterations, num_stages);

  for (auto* partition : subsystems)
  {
    cosim->addSubsystems(partition);
  }

  cosim->allocate();
  cosim->initializeSimulation(Tspan, &system->y(), &system->yp());

  const std::string filename = "cosim_solution_stages_"
                               + std::to_string(num_stages)
                               + "_iterations_" + std::to_string(max_iterations)
                               + "_intervals_" + std::to_string(num_intervals_1)
                               + "_" + std::to_string(num_intervals_2) + ".csv";

  std::ofstream output(filename);

  output << std::setprecision(16);

  auto out_cb = [&](const auto time, const auto& state) -> void
  {
    const auto* state_data = state.getData(GridKit::memory::HOST);

    output << time;

    for (size_t i = 0; i < static_cast<size_t>(state.getSize()) && i < MAX_STATE_OUTPUT; ++i)
    {
      output << "," << state_data[i];
    }

    output << '\n';
  };

  out_cb(0.0, system->y());
  // Run co-simulation and write the solution at every macro step.

  const auto t_start = std::chrono::steady_clock::now();

  cosim->runSimulation(out_cb);

  const auto   t_end     = std::chrono::steady_clock::now();
  const double elapsed_s = std::chrono::duration<double>(t_end - t_start).count();

  output << "\n# Metadata\n";
  output << "# solver: MultiStage Cosim\n";
  output << "# num_stages: " << num_stages << '\n';
  output << "# max_iterations: " << max_iterations << '\n';
  output << "# num_intervals_1: " << num_intervals_1 << '\n';
  output << "# num_intervals_2: " << num_intervals_2 << '\n';
  output << "# num_macro_steps: " << Tspan.size() - 1 << '\n';
  output << "# runtime_seconds: " << elapsed_s << '\n';
  output << "# num_partitions: " << num_partitions << '\n';

  output.close();

  std::cout << "\nCosimulation runtime: " << elapsed_s << " s  ("
            << Tspan.size() - 1 << " macro steps, "
            << 1e3 * elapsed_s / static_cast<double>(Tspan.size() - 1)
            << " ms/step)\n";

  // ---------------------------------------------------------------------------
  // Clean up
  // ---------------------------------------------------------------------------

  delete system;
  delete cosim;
  for (auto* partition : subsystems)
  {
    delete partition;
  }

  return 0;
}

template <typename ScalarT>
std::vector<ScalarT> linspace(ScalarT start, ScalarT end, std::size_t num_points)
{
  if (num_points == 0)
  {
    return {};
  }

  if (num_points == 1)
  {
    return {start};
  }

  std::vector<ScalarT> points(num_points);

  const ScalarT step = (end - start) / static_cast<ScalarT>(num_points - 1);

  for (std::size_t i = 0; i < num_points; ++i)
  {
    points[i] = start + static_cast<ScalarT>(i) * step;
  }

  // Guarantee the endpoint despite floating-point roundoff.
  points.back() = end;

  return points;
}