// SystemAssembler.hpp

#pragma once

#include <cassert>
#include <cmath>

#include <examples/PowerElectronics/ExamplesHelper/MicrogridNetwork.hpp>

/**
 * @brief Assemble a scaled microgrid network into a power electronics model.
 *
 * Adds the signal node, physical buses, generators, transmission lines,
 * loads, and virtual DQ buses stored in @p network to @p sys_model.
 *
 * This function does not construct or allocate any network components.
 *
 * @param[in] network Constructed scaled microgrid network whose components
 *                    are added to the system model.
 * @param[in,out] sys_model Power electronics model to which the network
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
 *       does not call PowerElectronicsModel::allocate().
 */
template <class ScalarT, typename IdxT>
void assembleSystem(ScaleMicrogridNetwork<ScalarT, IdxT>& network, GridKit::PowerElectronicsModel<ScalarT, IdxT>& sys_model)
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
