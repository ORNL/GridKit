#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
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

  constexpr int N_size  = 2;
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

  // --------------------------------------------------------------------------------
  // Create 2 Partitions and  a Partition interfaces to split to network
  // --------------------------------------------------------------------------------

  auto* bus_interface = new GridKit::PowerElectronics::BusPartitionInterface<double, int>(
      &network.buses[1],
      network.lines[1],
      14);

  if (int err = bus_interface->allocate())
  {
    return err;
  }

  auto* partition1 = new Subsystem();
  auto* partition2 = new Subsystem();

  // --------------------------------------------------------------------------------
  // Manually add components, nodes and a bus partition interface to Partition 1
  // --------------------------------------------------------------------------------

  partition1->addNode(&network.dg_signal);
  partition1->addComponent(network.generators[0]);
  partition1->addComponent(network.busesDQ[0]);
  partition1->addComponent(network.loads[0]);
  partition1->addNode(&network.buses[0]);
  partition1->addComponent(network.lines[0]);
  partition1->addComponent(network.generators[1]);
  partition1->addComponent(network.busesDQ[1]);
  partition1->addInterface(bus_interface);
  partition1->addNode(&network.buses[1]);

  // ---------------------------------------------------------------------------
  // Manually add components and nodes to Partition 2
  // ---------------------------------------------------------------------------

  partition2->addComponent(network.generators[2]);
  partition2->addComponent(network.busesDQ[2]);
  partition2->addComponent(network.loads[2]);
  partition2->addComponent(network.lines[1]);
  partition2->addNode(&network.buses[2]);
  partition2->addComponent(network.generators[3]);
  partition2->addComponent(network.busesDQ[3]);
  partition2->addComponent(network.lines[2]);
  partition2->addNode(&network.buses[3]);

  // ---------------------------------------------------------------------------
  // cosim
  // ---------------------------------------------------------------------------

  std::vector<double> time_points = linspace(0.0, 1.0, 200);
  // const double        dt          = 0.001;

  // std::vector<double> time_points = {
  //     0.0,
  //     dt,
  //     2.0 * dt};
  auto* cosim = new AnalysisManager::MultiStageCosimulation<double, int>();

  cosim->addSubsystems(partition1);
  cosim->addSubsystems(partition2);
  cosim->allocate();
  cosim->initializeSimulation(time_points, &system->y(), &system->yp());
  const auto t_start = std::chrono::steady_clock::now();

  cosim->runSimulation();

  const auto t_end = std::chrono::steady_clock::now();

  const double elapsed_s = std::chrono::duration<double>(t_end - t_start).count();

  std::cout << "\nCosimulation runtime: " << elapsed_s << " s  ("
            << time_points.size() - 1 << " macro steps, "
            << 1e3 * elapsed_s / static_cast<double>(time_points.size() - 1)
            << " ms/step)\n";

  auto* y_data = cosim->y().getData(GridKit::memory::HOST);

  std::cout << "\nFinal solution: y = [";

  for (int i = 0; i < cosim->y().getSize(); ++i)
  {
    std::cout << y_data[i];

    if (i + 1 < cosim->y().getSize())
    {
      std::cout << "\n";
    }
  }

  std::cout << "]\n";

  // ---------------------------------------------------------------------------
  // Clean up
  // ---------------------------------------------------------------------------
  delete system;
  delete partition1;
  delete partition2;

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
