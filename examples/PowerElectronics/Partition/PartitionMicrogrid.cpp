#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>

#include <GridKit/Model/PowerElectronics/PartitionInterface/BusLinePartitionInterface.hpp>
#include <GridKit/Model/PowerElectronics/SubsystemModel.hpp>
#include <GridKit/Model/PowerElectronics/SystemModel.hpp>
#include <GridKit/Solver/Dynamic/DynamicSolver.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>
#include <GridKit/Testing/Testing.hpp>

#include <examples/PowerElectronics/ExamplesHelper/MicrogridNetwork.hpp>
#include <examples/PowerElectronics/ExamplesHelper/PartitionUtilities.hpp>
#include <examples/PowerElectronics/ExamplesHelper/SystemAssembler.hpp>

using Component = GridKit::PowerElectronics::Component<double, size_t>;
using Node      = GridKit::PowerElectronics::NodeBase<double, size_t>;
using Subsystem = GridKit::PowerElectronics::SubsystemModel<double, size_t>;
using System    = GridKit::PowerElectronics::SystemModel<double, size_t>;

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

  constexpr size_t N_size = 2;
  constexpr double time   = 0.1;
  constexpr double alpha  = 0.1;

  bool use_jac = true;

  // ---------------------------------------------------------------------------
  // build the grid network and assemble the system model
  // ---------------------------------------------------------------------------
  ScaleMicrogridNetwork<double, size_t> network(N_size);

  auto* system = new System(use_jac);

  assembleSystemLeftToRight(network, *system);
  system->allocate();

  std::vector<double> y(system->size());
  std::vector<double> yp(system->size());

  for (size_t i = 0; i < system->size(); ++i)
  {
    y[i]  = static_cast<double>(i + 1);
    yp[i] = static_cast<double>(i + 1);
  }

  auto* system_y  = system->y().getData();
  auto* system_yp = system->yp().getData();

  for (size_t i = 0; i < system->size(); ++i)
  {
    system_y[i]  = y[i];
    system_yp[i] = yp[i];
  }

  system->y().setDataUpdated();
  system->yp().setDataUpdated();

  system->updateTime(time, alpha);
  system->evaluateResidual();
  system->evaluateJacobian();

  auto* system_jacobian = system->getCsrJacobian();
  auto* system_residual = system->getResidual().getData();

  // --------------------------------------------------------------------------------
  // Create 2 Partitions and  a Partition interfaces to split to network
  // --------------------------------------------------------------------------------

  auto* bus_interface = new GridKit::PowerElectronics::BusLinePartitionInterface<double, size_t>(
      &network.buses[1],
      network.lines[1],
      14);

  if (int err = bus_interface->allocate())
  {
    return err;
  }

  // --------------------------------------------------------------------------------
  // Manually add components, nodes and a bus partition interface to Partition 1
  // --------------------------------------------------------------------------------

  auto* partition1 = new Subsystem();

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

  auto* partition2 = new Subsystem();

  partition2->addComponent(network.generators[2]);
  partition2->addComponent(network.busesDQ[2]);
  partition2->addComponent(network.loads[2]);
  partition2->addComponent(network.lines[1]);
  partition2->addNode(&network.buses[2]);
  partition2->addComponent(network.generators[3]);
  partition2->addComponent(network.busesDQ[3]);
  partition2->addComponent(network.lines[2]);
  partition2->addNode(&network.buses[3]);

  partition1->allocate();
  partition2->allocate();

  std::vector<Subsystem*> partitions = {partition1, partition2};
  std::vector<double>     partition_residual(system->size(), 0.0);

  parallelFunctionEval(partitions, y, yp, partition_residual, time, alpha);

  // ---------------------------------------------------------------------------
  // Verify the subsystem Jacobians
  // ---------------------------------------------------------------------------

  double tolerance = 100 * std::numeric_limits<double>::epsilon();

  for (auto* partition : partitions)
  {
    partition->evaluateJacobian();

    bool jacobians_match = verifySubsystemJacobian(system_jacobian, *partition, tolerance);

    if (!jacobians_match)
    {
      std::cout << "ERROR: At least one subsystem Jacobian is incorrect!\n";
      return 1;
    }
  }

  // ---------------------------------------------------------------------------
  // Gather and verify the subsystem residuals against monolithic residuals
  // ---------------------------------------------------------------------------

  for (size_t i = 0; i < system->size(); ++i)
  {
    const bool matched = GridKit::Testing::isEqual(system_residual[i], partition_residual[i], tolerance);

    if (!matched)
    {
      std::cout << "ERROR: Residuals do not meet required tolerance!\n";
      return 1;
    }
  }

  // ---------------------------------------------------------------------------
  // Clean up
  // ---------------------------------------------------------------------------

  for (auto* partition : partitions)
  {
    delete partition;
  }

  delete system;

  return 0;
}
