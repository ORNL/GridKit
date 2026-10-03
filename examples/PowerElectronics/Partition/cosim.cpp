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
  yp[2]      = parms1.Vn_;
  yp[4]      = parms1.Kpv_ * parms1.Vn_;
  yp[6]      = (parms1.Kpc_ * parms1.Kpv_ * parms1.Vn_) / parms1.Lf_;
  yp[12 + 2] = parms1.Vn_;
  yp[12 + 4] = parms1.Kpv_ * parms1.Vn_;
  yp[12 + 6] = (parms1.Kpc_ * parms1.Kpv_ * parms1.Vn_) / parms1.Lf_;
  for (size_t i = 2; i < 4; i++)
  {
    yp[13 * i - 1 + 2] = parms2.Vn_;
    yp[13 * i - 1 + 4] = parms2.Kpv_ * parms2.Vn_;
    yp[13 * i - 1 + 6] = (parms2.Kpc_ * parms2.Kpv_ * parms2.Vn_) / parms2.Lf_;
  }

  // since the intial P_com = 0
  y[network.dg_signal.getNodeConnection(0).idx_] = parms1.wb_;

  system->y().setDataUpdated();
  system->yp().setDataUpdated();

  system->initialize();
  system->evaluateResidual();

  system->updateTime(0.0, 1.0e-8);
  system->evaluateJacobian();
  system->tagDifferentiable();

  auto* ida_ic = new AnalysisManager::Sundials::Ida<double, int>(system);

  ida_ic->setTolerance(1.0e-8);
  ida_ic->setMaxSteps(5000);
  ida_ic->configureSimulation();
  ida_ic->getDefaultInitialCondition();
  ida_ic->initializeSimulation(0.0, true);
  // ida_ic->runSimulation(0.0001);

  delete ida_ic;

  // --------------------------------------------------------------------------------
  // Create 2 Partitions and  a Partition interfaces to split to network
  // --------------------------------------------------------------------------------
  std::vector<Subsystem*> subsystems;
  partitionNetwork(network, subsystems, 48);

  // ---------------------------------------------------------------------------
  // cosim
  // ---------------------------------------------------------------------------

  std::vector<double> time_points  = linspace(0.0, 0.01, 10);
  std::vector<double> time_points1 = linspace(0.01, 1.0, 10);

  time_points.insert(time_points.end(), time_points1.begin() + 1, time_points1.end());

  // std::vector<double> time_points_b = linspace(1.0, 1.01, 50);
  // std::vector<double> time_points2  = linspace(1.01, 2.0, 40);
  // time_points_b.insert(time_points_b.end(), time_points2.begin() + 1, time_points2.end());

  auto* cosim = new AnalysisManager::MultiStageCosimulation<double, int>(2, 4);

  for (auto* partition : subsystems)
  {
    cosim->addSubsystems(partition);
  }

  cosim->allocate();
  cosim->initializeSimulation(time_points, &system->y(), &system->yp());
  const auto t_start = std::chrono::steady_clock::now();

  std::ofstream output("cosim_solution.csv");

  output << std::setprecision(16);

  // Header
  output << "time";

  for (int i = 0; i < cosim->y().getSize() && i < 500; ++i)
  {
    output << ",y" << i;
  }

  output << '\n';

  // Initial condition
  const auto* y_initial = cosim->y().getData(GridKit::memory::HOST);

  output << time_points.front();

  for (int i = 0; i < cosim->y().getSize() && i < 500; ++i)
  {
    output << ',' << y_initial[i];
  }

  output << '\n';

  auto out_cb = [&](const auto time, const auto& state) -> void
  {
    const auto* state_data = state.getData(GridKit::memory::HOST);

    output << time;

    for (int i = 0; i < state.getSize() && i < 500; ++i)
    {
      output << ',' << state_data[i];
    }

    output << '\n';
  };

  // Run co-simulation and write the solution at every macro step.
  cosim->runSimulation(out_cb);

  // network.loads[0]->setResistance(3.0);
  // network.loads[0]->setInductance(6.4e-3);

  // network.loads[2]->setResistance(2.0);
  // network.loads[2]->setInductance(3.2e-3);

  // cosim->initializeSimulation(time_points_b, &cosim->y(), &cosim->yp());
  // cosim->runSimulation(out_cb);

  output.close();

  const auto t_end = std::chrono::steady_clock::now();

  const double elapsed_s = std::chrono::duration<double>(t_end - t_start).count();

  std::cout << "\nCosimulation runtime: " << elapsed_s << " s  ("
            << time_points.size() - 1 << " macro steps, "
            << 1e3 * elapsed_s / static_cast<double>(time_points.size() - 1)
            << " ms/step)\n";

  // auto* y_data = cosim->y().getData(GridKit::memory::HOST);

  // std::cout << "\nFinal solution: y = [";

  // for (int i = 0; i < cosim->y().getSize(); ++i)
  // {
  //   std::cout << y_data[i];

  //   if (i + 1 < cosim->y().getSize())
  //   {
  //     std::cout << "\n";
  //   }
  // }

  // std::cout << "]\n";

  // ---------------------------------------------------------------------------
  // Clean up
  // ---------------------------------------------------------------------------
  delete system;
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
