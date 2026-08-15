#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <vector>

#include <GridKit/Model/PowerElectronics/Bus/MicrogridBus.hpp>
#include <GridKit/Model/PowerElectronics/Bus/SignalNode.hpp>
#include <GridKit/Model/PowerElectronics/DistributedGenerator/DistributedGenerator.hpp>
#include <GridKit/Model/PowerElectronics/MicrogridBusDQ/MicrogridBusDQ.hpp>
#include <GridKit/Model/PowerElectronics/MicrogridLine/MicrogridLine.hpp>
#include <GridKit/Model/PowerElectronics/MicrogridLoad/MicrogridLoad.hpp>
#include <GridKit/Model/PowerElectronics/PartitionInterface/BusPartitionInterface.hpp>
#include <GridKit/Model/PowerElectronics/SubsystemModel.hpp>
#include <GridKit/Model/PowerElectronics/SystemModelPowerElectronics.hpp>
#include <GridKit/Solver/Dynamic/DynamicSolver.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>

#include "Common/JacTestHelper.hpp"
#include "Common/MicrogridNetwork.hpp"
#include "Common/PartitionUtilities.hpp"

using Component = GridKit::CircuitComponent<double, size_t>;
using Node      = GridKit::PowerElectronics::NodeBase<double, size_t>;
using Subsystem = GridKit::SubsystemModel<double, size_t>;
using System    = GridKit::PowerElectronicsModel<double, size_t>;

std::vector<size_t> getComponentConnections(const std::vector<Component*>& components);
std::vector<size_t> getNodeConnections(const std::vector<Node*>& nodes);

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

  constexpr size_t num_network_sections = 2;
  constexpr double time                 = 0.1;
  constexpr double alpha                = 0.1;

  bool use_jac = true;

  // ---------------------------------------------------------------------------
  // build the grid network and assemble the system model
  // ---------------------------------------------------------------------------
  GridKit::ScaleMicrogridNetwork network(num_network_sections);

  GridKit::buildScaleMicrogridNetwork(network);

  auto* system = new System(use_jac);

  GridKit::assembleSystemLeftToRight(network, *system);
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

  //------------------------------------------------------------------------------
  // Gather all global indices belonging partition 1 and 2 to test release() later
  //------------------------------------------------------------------------------
  std::vector<Component*> components = {
      network.generators[0],
      network.generators[1],
      network.lines[1],
      network.loads[0],
      network.busesDQ[0],
      network.busesDQ[1],
      network.generators[2],
      network.generators[3],
      network.lines[2],
      network.lines[3],
      network.loads[2],
      network.busesDQ[2],
      network.busesDQ[3]};

  std::vector<Node*> nodes = {
      &network.dg_signal,
      &network.buses[0],
      &network.buses[1],
      &network.buses[2],
      &network.buses[3]};

  const auto original_component_connections = getComponentConnections(components);
  const auto original_node_connections      = getNodeConnections(nodes);

  // --------------------------------------------------------------------------------
  // Create 2 Partitions and Partition interfaces
  // --------------------------------------------------------------------------------

  auto* bus_interface = new GridKit::BusPartitionInterface<double, size_t>(
      &network.buses[1],
      network.lines[2],
      14);

  bus_interface->allocate();

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
  partition1->addComponent(network.lines[1]);
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
  partition2->addComponent(network.lines[2]);
  partition2->addNode(&network.buses[2]);
  partition2->addComponent(network.generators[3]);
  partition2->addComponent(network.busesDQ[3]);
  partition2->addComponent(network.lines[3]);
  partition2->addNode(&network.buses[3]);

  std::vector<Subsystem*> partitions = {partition1, partition2};

  for (auto* partition : partitions)
  {
    partition->allocate();
  }

  std::vector<double> partition_residual(system->size(), 0.0);

  GridKit::evaluatePartitionResiduals(partitions, y, yp, partition_residual, time, alpha);

  // ---------------------------------------------------------------------------
  // Verify the subsystem Jacobians
  // ---------------------------------------------------------------------------

  bool jacobians_match = true;

  for (auto* partition : partitions)
  {
    partition->evaluateJacobian();

    auto* partition_jacobian = partition->getCsrJacobian();

    jacobians_match = jacobians_match && GridKit::Testing::verifySubsystemJacobian(*system_jacobian, *partition_jacobian, *partition);
  }

  if (!jacobians_match)
  {
    std::cout << "ERROR: At least one subsystem Jacobian is incorrect!\n";
    return 1;
  }

  // ---------------------------------------------------------------------------
  // Gather and verify the subsystem residuals against monolithic residuals
  // ---------------------------------------------------------------------------

  double max_error = 0.0;

  for (size_t i = 0; i < system->size(); ++i)
  {
    const double error = std::abs(system_residual[i] - partition_residual[i]) / std::abs(system_residual[i] + 1);

    max_error = std::max(max_error, error);
  }

  const double machine_epsilon =
      std::numeric_limits<double>::epsilon();

  const bool residuals_match = max_error <= machine_epsilon;

  std::cout << "\nPartition Microgrid Validation\n";
  std::cout << "------------------------------\n";

  std::cout << std::left
            << std::setw(32) << "Maximum residual error:"
            << std::setprecision(16)
            << max_error << '\n';

  std::cout << std::left
            << std::setw(32) << "Machine epsilon:"
            << machine_epsilon << '\n';

  std::cout << std::left
            << std::setw(32) << "Residuals matched:"
            << (residuals_match ? "True" : "False")
            << '\n';

  // ---------------------------------------------------------------------------
  // Verify subsystem release() from SubsystemModel
  // ---------------------------------------------------------------------------

  for (auto* partition : partitions)
  {
    partition->release();
  }

  const bool components_restored = getComponentConnections(components) == original_component_connections;
  const bool nodes_restored      = getNodeConnections(nodes) == original_node_connections;

  if (!components_restored || !nodes_restored)
  {
    std::cout << "ERROR: Subsystem release did not restore "
                 "the original global connection indices!\n";

    return 1;
  }

  // ---------------------------------------------------------------------------
  // Clean up
  // ---------------------------------------------------------------------------
  delete system;

  for (auto* partition : partitions)
  {
    delete partition;
  }

  return residuals_match ? 0 : 1;
}

/**
 * @brief Collect the connection indices of a set of components.
 */
std::vector<size_t> getComponentConnections(const std::vector<Component*>& components)
{
  std::vector<size_t> connections;

  for (const auto* component : components)
  {
    for (size_t i = 0; i < component->size(); ++i)
    {
      connections.push_back(component->getNodeConnection(i));
    }
  }

  return connections;
}

/**
 * @brief Collect the connection indices of a set of nodes.
 */
std::vector<size_t> getNodeConnections(const std::vector<Node*>& nodes)
{
  std::vector<size_t> connections;

  for (auto* node : nodes)

    for (size_t i = 0; i < node->size(); ++i)
    {
      connections.push_back(node->getNodeConnection(i).idx_);
    }

  return connections;
}