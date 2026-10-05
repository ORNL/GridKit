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


struct GenDataIndex {
  bool refframe;
  int idx_y8;
  int idx_y9;
  int idx_y10;
  int idx_y11;
  int idx_y12;
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

  constexpr int N_size  = 5000;
  bool          use_jac = true;

  // ---------------------------------------------------------------------------
  // build the grid network and assemble the system model
  // ---------------------------------------------------------------------------
  ScaleMicrogridNetwork<double, int> network(N_size);

  auto* system = new System(use_jac);

  assembleSystem(network, *system);
  system->allocate();

  std::vector<GenDataIndex> genData(5);

//Identify the location for the variables need to compute P and Q for each generator
  for(int i=0; i < 5; i++)
  {
    auto* gen  = network.generators[i];
    auto& gendata = genData[i];

    gendata.refframe = (i == 0);
    gendata.idx_y8   = gen->getNodeConnection(8);
    gendata.idx_y9   = gen->getNodeConnection(9);
    gendata.idx_y10  = gen->getNodeConnection(10);
    gendata.idx_y11  = gen->getNodeConnection(11);
    gendata.idx_y12  = gendata.refframe? -1:gen->getNodeConnection(12);
  }

  // Generator parameters used to construct the initial conditions.
  const auto& parms1 = network.DGParam_list[0];
  const auto& parms2 = network.DGParam_list[2];

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
  // DGs 1
  for (int i = 0; i < 2 * N_size; ++i)
  {
    const auto& params = network.DGParam_list[i];

    yp[13 * i - 1 + 2] = params.Vn_;
    yp[13 * i - 1 + 4] = params.Kpv_ * params.Vn_;
    yp[13 * i - 1 + 6] = (params.Kpc_ * params.Kpv_ * params.Vn_) / params.Lf_;
  }

  // Initial P_com = 0, so initialize the reference-frame frequency.
  y[network.dg_signal.getNodeConnection(0).idx_] = network.DGParam_list[0].wb_;

  system->y().setDataUpdated();
  system->yp().setDataUpdated();

  // --------------------------------------------------------------------------------
  // Create 2 Partitions and  a Partition interfaces to split to network
  // --------------------------------------------------------------------------------
  std::vector<Subsystem*> subsystems;
  partitionNetwork(network, subsystems, 1000);

  // ---------------------------------------------------------------------------
  // cosim
  // ---------------------------------------------------------------------------
  constexpr std::size_t max_iterations = 1;
  constexpr std::size_t num_stages     = 1;

  constexpr std::size_t num_intervals_1 = 140;
  constexpr std::size_t num_intervals_2 = 0;
  // Intervals to use
  // 15 30
  // 20 40
  // 25 50
  // 30 60
  // 35 70

  std::vector<double> time_points  = linspace(0.0, 1.0, num_intervals_1 + 1);
  // std::vector<double> time_points1 = linspace(0.01, 1.0, num_intervals_2 + 1);

  // time_points.insert(time_points.end(), time_points1.begin() + 1, time_points1.end());

  auto* cosim = new AnalysisManager::MultiStageCosimulation<double, int>(max_iterations, num_stages);

  for (auto* partition : subsystems)
  {
    cosim->addSubsystems(partition);
  }

  cosim->allocate();
  cosim->initializeSimulation(time_points, &system->y(), &system->yp());

  const std::string filename = "cosim_solution_stages_" + std::to_string(num_stages) +
                               "_iterations_" + std::to_string(max_iterations) +
                               "_intervals_" + std::to_string(num_intervals_1) +
                               "_" + std::to_string(num_intervals_2) + ".csv";

  std::ofstream output(filename);

  output << std::setprecision(16);

  auto out_cb = [&](const auto time, const auto& state) -> void
  {
    const auto* state_data = state.getData(GridKit::memory::HOST);

    output << time;

    for (int i = 0; i < static_cast<int>(genData.size()); ++i)
    {
    const auto& gendata = genData[i];

    const auto y8  = state_data[gendata.idx_y8];
    const auto y9  = state_data[gendata.idx_y9];
    const auto y10 = state_data[gendata.idx_y10];
    const auto y11 = state_data[gendata.idx_y11];

    const auto delta = gendata.refframe ? 0.0 : state_data[gendata.idx_y12];
    const auto P     = y8 * y10 + y9 * y11;
    const auto Q     = -y8 * y11 + y9 * y10;

    output << ";   " << delta
           << ", "   << P
           << ", "   << Q;
    }

    output << '\n';
  };

  out_cb(0.0, system->y());
  // Run co-simulation and write the solution at every macro step.
  

  const auto t_start = std::chrono::steady_clock::now(); 

  cosim->runSimulation(out_cb);

  const auto t_end = std::chrono::steady_clock::now();
  const double elapsed_s = std::chrono::duration<double>(t_end - t_start).count();

  output << "\n# Metadata\n";
  output << "# solver: MultiStage Cosim\n";
  output << "# num_stages: "      << num_stages << '\n';
  output << "# max_iterations: "  << max_iterations << '\n';
  output << "# num_intervals_1: " << num_intervals_1 << '\n';
  output << "# num_intervals_2: " << num_intervals_2 << '\n';
  output << "# num_macro_steps: " << time_points.size() - 1 << '\n';
  output << "# runtime_seconds: " << elapsed_s << '\n';

  output.close();

  std::cout << "\nCosimulation runtime: " << elapsed_s << " s  ("
            << time_points.size() - 1 << " macro steps, "
            << 1e3 * elapsed_s / static_cast<double>(time_points.size() - 1)
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