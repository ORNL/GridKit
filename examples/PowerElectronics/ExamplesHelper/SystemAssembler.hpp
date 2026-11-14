// SystemAssembler.hpp

#pragma once

#include <cassert>
#include <cmath>

#include <examples/PowerElectronics/ExamplesHelper/MicrogridNetwork.hpp>

/**
 * @brief Adds Microgrid Network components into system model.
 *
 * Adds the signal node, physical buses, generators, transmission lines,
 * loads, and virtual DQ buses stored in @p network to @p sys_model.
 *
 * @param[in] network Constructed scaled microgrid network whose components
 *                    are added to the system model.
 * @param[in,out] sys_model system model to which the network
 *                          components and nodes are added.
 *
 * @pre @c network.N_size is greater than zero.
 * @pre All component and node pointers referenced by @p network are valid.
 *
 * @post The signal node and all physical buses in @p network have been added
 *       to @p sys_model.
 * @post All generators, transmission lines, loads, and virtual DQ buses in
 *       @p network have been added to @p sys_model.
 *
 * @note This function only assembles the network into the system model. It
 *       does not call SystemModel::allocate().
 */
template <class ScalarT, typename IdxT>
void assembleSystem(ScaleMicrogridNetwork<ScalarT, IdxT>&                  network,
                    GridKit::PowerElectronics::SystemModel<ScalarT, IdxT>& sys_model)
{
  size_t N_size = network.N_size;

  // Ensure minimum size requirement
  assert(N_size > 0);

  // Add all bus nodes
  sys_model.addNode(&network.dg_signal);

  for (size_t i = 0; i < 2 * N_size; i++)
  {
    sys_model.addNode(&network.buses[i]);
  }

  // Add all generators
  for (IdxT i = 0; i < 2 * N_size; i++)
  {
    sys_model.addComponent(network.generators[i]);
  }

  // Load all the Line components
  for (IdxT i = 0; i < 2 * N_size - 1; i++)
  {
    sys_model.addComponent(network.lines[i]);
  }

  //  Load all the Load components
  for (IdxT i = 0; i < 2 * N_size; i++)
  {
    if (network.loads[i] != nullptr)
    {
      sys_model.addComponent(network.loads[i]);
    }
  }

  // Add all the microgrid Virtual DQ Buses
  for (IdxT i = 0; i < 2 * N_size; i++)
  {
    sys_model.addComponent(network.busesDQ[i]);
  }
}

/**
 * @brief Assemble the scaled microgrid from left to right.
 *
 * Adds the network components to @p sys_model in physical left-to-right
 * order. For each bus location, the generator, virtual DQ bus, load, line,
 * and bus associated with that location are added before moving to the next
 * location in the network.
 *
 * This ordering is useful when comparing the monolithic system with
 * partition-based evaluations that traverse the network from left to right.
 *
 * @param[in] network Constructed scaled microgrid network.
 * @param[in,out] sys_model Power electronics model to populate.
 *
 * @pre @c network.N_size is greater than zero.
 * @pre @p network has been constructed by buildScaleMicrogridNetwork().
 * @pre All component and node pointers stored in @p network are valid.
 *
 * @post All nodes and components in @p network have been added to
 *       @p sys_model in left-to-right network order.
 */
template <class ScalarT, typename IdxT>
void assembleSystemLeftToRight(
    ScaleMicrogridNetwork<ScalarT, IdxT>&                  network,
    GridKit::PowerElectronics::SystemModel<ScalarT, IdxT>& sys_model)
{
  const size_t num_ibrs = 2 * network.N_size;

  assert(network.N_size > 0);

  sys_model.addNode(&network.dg_signal);

  for (IdxT i = 0; i < num_ibrs; ++i)
  {
    sys_model.addComponent(network.generators[i]);
    sys_model.addComponent(network.busesDQ[i]);

    if (network.loads[i] != nullptr)
    {
      sys_model.addComponent(network.loads[i]);
    }

    if (i > 0)
    {
      sys_model.addComponent(network.lines[i - 1]);
    }

    sys_model.addNode(&network.buses[i]);
  }
}
