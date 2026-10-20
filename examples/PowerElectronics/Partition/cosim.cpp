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

struct GenDataIndex
{
  bool refframe;
  int  idx_y11;
  int  idx_y12;
  int  idx_y13;
  int  idx_y14;
  int  idx_y15;
};

template <typename ScalarT>
std::vector<ScalarT> linspace(ScalarT start, ScalarT end, std::size_t num_points);
int runCosim(size_t num_intervals_1, 
             size_t num_intervals_2,
             size_t max_iterations,
             size_t num_stages,
             size_t num_partitions, 
             double tol, 
             size_t N_size, 
             std::string type);


/*
 * Verify partitioned residual and Jacobian evaluation against the
 *        monolithic microgrid model.
 *
 * The microgrid is manually divided into two subsystems. A partition interface
 * is introduced at the boundary between bus 1 and the line connecting buses 1
 * and 2. Each subsystem is evaluated independently, and its residuals are gathered
 * and compared with the monolithic reference.
 */

 constexpr bool WRITE_OUTPUT = true;

 // Those generators whose data we want to output
 const std::vector<size_t> generator_ids = {9, 10, 11, 99, 9999};

 int main()
 {
    std::string timing = "timing" ;

    if(WRITE_OUTPUT)
    { 
      timing = ""; //If output is to be written, we do not consider timing
    } 

    std::string info_type = "_tol_1e-6_" + timing;
   // MS1_1(ZOH)
    runCosim(5,  20, 1, 1, 100, 1e-6, 5000, info_type);
    runCosim(10, 30, 1, 1, 100, 1e-6, 5000, info_type);
    runCosim(15, 40, 1, 1, 100, 1e-6, 5000, info_type);
    runCosim(20, 50, 1, 1, 100, 1e-6, 5000, info_type);
    runCosim(25, 60, 1, 1, 100, 1e-6, 5000, info_type);
    runCosim(30, 70, 1, 1, 100, 1e-6, 5000, info_type);
    runCosim(35, 80, 1, 1, 100, 1e-6, 5000, info_type);

    // MS1_2
    runCosim(5,  20, 1, 2, 100, 1e-6, 5000, info_type);
    runCosim(10, 30, 1, 2, 100, 1e-6, 5000, info_type);
    runCosim(15, 40, 1, 2, 100, 1e-6, 5000, info_type);
    runCosim(20, 50, 1, 2, 100, 1e-6, 5000, info_type);
    runCosim(25, 60, 1, 2, 100, 1e-6, 5000, info_type);
    runCosim(30, 70, 1, 2, 100, 1e-6, 5000, info_type);
    runCosim(35, 80, 1, 2, 100, 1e-6, 5000, info_type);

    //MS2_2
    runCosim(5,  20, 2, 2, 100, 1e-6, 5000, info_type);
    runCosim(10, 30, 2, 2, 100, 1e-6, 5000, info_type);
    runCosim(15, 40, 2, 2, 100, 1e-6, 5000, info_type);
    runCosim(20, 50, 2, 2, 100, 1e-6, 5000, info_type);
    runCosim(25, 60, 2, 2, 100, 1e-6, 5000, info_type);
    runCosim(30, 70, 2, 2, 100, 1e-6, 5000, info_type);
    runCosim(35, 80, 2, 2, 100, 1e-6, 5000, info_type);

    // MS1_3
    runCosim(5,  20, 1, 3, 100, 1e-6, 5000, info_type);
    runCosim(10, 30, 1, 3, 100, 1e-6, 5000, info_type);
    runCosim(15, 40, 1, 3, 100, 1e-6, 5000, info_type);
    runCosim(20, 50, 1, 3, 100, 1e-6, 5000, info_type);
    runCosim(25, 60, 1, 3, 100, 1e-6, 5000, info_type);
    runCosim(30, 70, 1, 3, 100, 1e-6, 5000, info_type);
    runCosim(35, 80, 1, 3, 100, 1e-6, 5000, info_type);

    // MS2_3
    runCosim(5,  20, 2, 3, 100, 1e-6, 5000, info_type);
    runCosim(10, 30, 2, 3, 100, 1e-6, 5000, info_type);
    runCosim(15, 40, 2, 3, 100, 1e-6, 5000, info_type);
    runCosim(20, 50, 2, 3, 100, 1e-6, 5000, info_type);
    runCosim(25, 60, 2, 3, 100, 1e-6, 5000, info_type);
    runCosim(30, 70, 2, 3, 100, 1e-6, 5000, info_type);
    runCosim(35, 80, 2, 3, 100, 1e-6, 5000, info_type);

    // MS1_4
    runCosim(5,  20, 1, 4, 100, 1e-6, 5000, info_type);
    runCosim(10, 30, 1, 4, 100, 1e-6, 5000, info_type);
    runCosim(15, 40, 1, 4, 100, 1e-6, 5000, info_type);
    runCosim(20, 50, 1, 4, 100, 1e-6, 5000, info_type);
    runCosim(25, 60, 1, 4, 100, 1e-6, 5000, info_type);
    runCosim(30, 70, 1, 4, 100, 1e-6, 5000, info_type);
    runCosim(35, 80, 1, 4, 100, 1e-6, 5000, info_type);

    // MS2_4
    runCosim(5,  20, 2, 4, 100, 1e-6, 5000, info_type);
    runCosim(10, 30, 2, 4, 100, 1e-6, 5000, info_type);
    runCosim(15, 40, 2, 4, 100, 1e-6, 5000, info_type);
    runCosim(20, 50, 2, 4, 100, 1e-6, 5000, info_type);
    runCosim(25, 60, 2, 4, 100, 1e-6, 5000, info_type);
    runCosim(30, 70, 2, 4, 100, 1e-6, 5000, info_type);
    runCosim(35, 80, 2, 4, 100, 1e-6, 5000, info_type);

    return 0;

 }

int runCosim(size_t num_intervals_1, 
             size_t num_intervals_2,
             size_t max_iterations,
             size_t num_stages,
             size_t num_partitions, 
             double tol, 
             size_t N_size, 
             std::string type)
{

  constexpr size_t load_faulted_index = 10;
  bool             use_jac          = true;

  // ---------------------------------------------------------------------------
  // build the grid network and assemble the system model
  // ---------------------------------------------------------------------------
  ScaleMicrogridNetwork<double, int> network(N_size);

  auto* system = new System(use_jac);

  assembleSystem(network, *system);
  system->allocate();

  // ---------------------------------------------------------------------------
  // Collect the indices of the generators to inspect
  // ---------------------------------------------------------------------------
  std::vector<GenDataIndex> gen_data;

  for (const auto i : generator_ids)
  {
    const auto* gen = network.generators[i];
    const bool refframe = (i == 0);

    gen_data.push_back({
        refframe,
        gen->getNodeConnection(11),
        gen->getNodeConnection(12),
        gen->getNodeConnection(13),
        gen->getNodeConnection(14),
        refframe ? -1 : gen->getNodeConnection(15)
    });
  }

  // Apply a fault to the selected load
  for (size_t i = 0; i < network.loads.size(); i += 2)
  {
    if(i == load_faulted_index)
    {
      network.loads[i]->setResistance(3);
      network.loads[i]->setInductance(0.1);
    }
    else
    {
      network.loads[i]->setResistance(5);
      network.loads[i]->setInductance(20);
    }
  }

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

  // ---------------------------------------------------------------------------
  // Check initial condition consistency
  // ---------------------------------------------------------------------------
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
  // Create Partitions and  a Partition interfaces to split to network
  // --------------------------------------------------------------------------------
  std::vector<Subsystem*> subsystems;
  partitionNetwork(network, subsystems, num_partitions);

  // ---------------------------------------------------------------------------
  // Configure co-simulation
  // ---------------------------------------------------------------------------

  std::vector<double> Tspan1 = linspace(0.0, 0.01, num_intervals_1 + 1);
  std::vector<double> Tspan2 = linspace(0.01, 1.0, num_intervals_2 + 1);
  std::vector<double> Tspan  = Tspan1;

  Tspan.insert(Tspan.end(), Tspan2.begin() + 1, Tspan2.end());

  auto* cosim = new AnalysisManager::MultiStageCosimulation<double, int>(max_iterations, num_stages, tol);

  for (auto* partition : subsystems)
  {
    cosim->addSubsystems(partition);
  }

  cosim->allocate();
  cosim->initializeSimulation(Tspan, &system->y(), &system->yp());

  const std::string filename = "cosim_sol_stages_"
                               + std::to_string(num_stages)
                               + "_iterations_" + std::to_string(max_iterations)
                               + "_intervals_" + std::to_string(num_intervals_1)
                               + "_" + std::to_string(num_intervals_2)
                               + "_partitions_" + std::to_string(num_partitions) 
                               + type +  ".csv";

  std::ofstream output(filename);

  output << std::setprecision(16);

  if constexpr (WRITE_OUTPUT)
  {
    output << "time";

    for (const auto i : generator_ids)
    {
      output << ",delta_" << i
            << ",P_" << i
            << ",Q_" << i;
    }

    output << '\n';
  }

  auto out_cb = [&](const auto time, const auto& state) -> void
  {
    const auto* state_data = state.getData(GridKit::memory::HOST);

    output << time;

    for (const auto& gen : gen_data)
    {
      const auto y11 = state_data[gen.idx_y11];
      const auto y12 = state_data[gen.idx_y12];
      const auto y13 = state_data[gen.idx_y13];
      const auto y14 = state_data[gen.idx_y14];

      const auto delta = gen.refframe ? 0.0 : state_data[gen.idx_y15];
      const auto P     = y11 * y13 + y12 * y14;
      const auto Q     = -y11 * y14 + y12 * y13;

      output << "," << delta
            << "," << P
            << "," << Q;
    }

    output << '\n';
  };

  if constexpr (WRITE_OUTPUT)
  {
    out_cb(0.0, system->y());
  }
  // Run co-simulation and write the solution at every macro step.

  const auto t_start = std::chrono::steady_clock::now();

  if constexpr (WRITE_OUTPUT)
  {
    cosim->runSimulation(out_cb);
  }
  else
  {
    cosim->runSimulation();
  }

  const auto   t_end     = std::chrono::steady_clock::now();
  const double elapsed_s = std::chrono::duration<double>(t_end - t_start).count();

  if constexpr (!WRITE_OUTPUT)
  {
    output << "\n# Metadata\n";
    output << "# solver: MultiStage Cosim\n";
    output << "# num_stages: " << num_stages << '\n';
    output << "# max_iterations: " << max_iterations << '\n';
    output << "# num_intervals_1: " << num_intervals_1 << '\n';
    output << "# num_intervals_2: " << num_intervals_2 << '\n';
    output << "# num_macro_steps: " << Tspan.size() - 1 << '\n';
    output << "# num_partitions: " << num_partitions << '\n';
    output << "# runtime_seconds: " << elapsed_s  << '\n';
  }

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