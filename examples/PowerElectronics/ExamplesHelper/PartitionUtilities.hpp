// PartitionUtilities.hpp
#pragma once

#include <omp.h>

#include <cassert>
#include <cstddef>
#include <vector>

#include <GridKit/Model/PowerElectronics/PartitionInterface/BusPartitionInterface.hpp>
#include <GridKit/Model/PowerElectronics/SubsystemModel.hpp>

#include <examples/PowerElectronics/ExamplesHelper/MicrogridNetwork.hpp>

/**
 * @brief Partition a scaled microgrid network into subsystem models.
 *
 * Divides the network into contiguous groups of IBRs and creates one
 * subsystem for each group. The IBRs are distributed as evenly as possible,
 * with any remainder assigned to the first partitions.
 *
 * The reference signal node is added to the first subsystem. A partition
 * interface is added at the right boundary of every subsystem except the
 * last to represent its connection to the neighboring subsystem.
 *
 * @param[in,out] network Scaled microgrid network to partition.
 * @param[out] subsystems Vector populated with the created subsystem models.
 * @param[in] num_partitions Number of subsystems to create.
 *
 * @pre The network has been constructed and contains @c 2*N_size IBRs.
 * @pre @p num_partitions is greater than zero and does not exceed the
 *      number of IBRs in the network.
 * @pre The generators, buses, virtual DQ buses, loads, and lines referenced
 *      by the network have valid lifetimes for use by the created subsystems.
 *
 * @post @p subsystems contains exactly @p num_partitions subsystem models.
 * @post Every IBR in the network belongs to exactly one subsystem.
 * @post Every subsystem except the last has a partition interface at its
 *       right boundary.
 *
 * @note The subsystem models are dynamically allocated. The caller is
 *       responsible for releasing and deleting them.
 */
template <class ScalarT, typename IdxT>
void partitionNetwork(
    ScaleMicrogridNetwork<ScalarT, IdxT>&                 network,
    std::vector<GridKit::SubsystemModel<ScalarT, IdxT>*>& subsystems,
    size_t                                                num_partitions)
{
  const size_t num_ibrs = 2 * network.N_size;

  assert(num_partitions <= num_ibrs);

  subsystems.resize(num_partitions);

  IdxT q = num_ibrs / num_partitions;
  IdxT r = num_ibrs % num_partitions;

  IdxT partition_begin = 0;

  for (IdxT j = 0; j < num_partitions; j++)
  {
    auto* partition = new GridKit::SubsystemModel<ScalarT, IdxT>();

    // Add the reference signal node to the first partition.
    if (j == 0)
    {
      partition->addNode(&network.dg_signal);
    }

    // Divide the IBRs as evenly as possible among the partitions.
    // Each partition receives q IBRs, and the first r partitions receive
    // one additional IBR to account for any remainder. The partition spans
    // [partition_begin, partition_end).
    const IdxT partition_size = q + ((j < r) ? 1 : 0);
    const IdxT partition_end  = partition_begin + partition_size;

    // Add all components belonging to this partition.
    for (IdxT i = partition_begin; i < partition_end; ++i)
    {
      partition->addComponent(network.generators[i]);
      partition->addComponent(network.busesDQ[i]);

      if (network.loads[i] != nullptr)
      {
        partition->addComponent(network.loads[i]);
      }

      // Each line is owned by the partition containing its right bus.
      // line[i - 1] connects bus[i - 1] to bus[i].
      if (i > 0)
      {
        partition->addComponent(network.lines[i - 1]);
      }

      partition->addNode(&network.buses[i]);
    }

    // If another partition exists to the right, the line crossing this
    // partition's right boundary is owned by that next partition. Add an
    // interface here to the left partition.
    if (partition_end < num_ibrs)
    {

      const IdxT boundary_line_index = partition_end - 1;

      auto* bus_interface = new GridKit::BusPartitionInterface<ScalarT, IdxT>(
          &network.buses[boundary_line_index],
          network.lines[boundary_line_index],
          network.model_id_next++);

      bus_interface->allocate();
      partition->addInterface(bus_interface);
    }

    subsystems[j]   = partition;
    partition_begin = partition_end;
  }
}

/**
 * @brief Evaluate subsystem residuals and reconstruct the global residual.
 *
 * Distributes the global state and state-derivative vectors to each
 * subsystem using its internal and external index mappings. Each subsystem
 * residual is then evaluated in parallel, and its internal residual entries
 * are gathered into the global residual vector.
 *
 * @param[in] subsystems Subsystem models to evaluate.
 * @param[in] y Global state vector.
 * @param[in] yp Global state-derivative vector.
 * @param[out] f Global residual vector reconstructed from the subsystem
 *               residuals.
 * @param[in] time Current simulation time.
 * @param[in] alpha Jacobian scaling parameter associated with the current
 *                  time-integration evaluation.
 *
 * @pre All subsystem models have been allocated.
 * @pre @p y and @p yp contain all global entries
 * @pre @p f is large enough to contain every global residual entry referenced
 *      by the subsystems.
 * @pre The subsystems have disjoint internal index sets so that parallel
 *      writes to @p f do not overlap.
 *
 * @post Each subsystem residual has been evaluated at the supplied
 *       @p time and @p alpha.
 * @post @p f contains the reconstructed residual for all internal variables
 *       represented by the subsystems.
 */
template <class ScalarT, typename IdxT>
void evaluatePartitionResiduals(
    const std::vector<GridKit::SubsystemModel<ScalarT, IdxT>*>& subsystems,
    const std::vector<ScalarT>&                                 y,
    const std::vector<ScalarT>&                                 yp,
    std::vector<ScalarT>&                                       f,
    ScalarT                                                     time,
    ScalarT                                                     alpha)
{
#ifdef _OPENMP
#pragma omp parallel for schedule(guided)
#endif
  for (auto* partition : subsystems)
  {
    partition->updateTime(time, alpha);

    auto* external_y  = partition->getExternalDataY().getData();
    auto* external_yp = partition->getExternalDataYP().getData();

    // Supply external variable values required by this partition from neighboring subsystems.
    for (size_t i = 0; i < partition->getExternSize(); ++i)
    {
      const auto global_index = partition->getExternalDataIndices()[i];

      external_y[i]  = y[global_index];
      external_yp[i] = yp[global_index];
    }

    partition->getExternalDataY().setDataUpdated();
    partition->getExternalDataYP().setDataUpdated();

    auto* partition_y  = partition->y().getData();
    auto* partition_yp = partition->yp().getData();

    // Supply external variable values required by this partition.
    for (size_t i = 0; i < partition->getInternalSize(); i++)
    {
      partition_y[i]  = y[partition->getNodeConnection(i)];
      partition_yp[i] = yp[partition->getNodeConnection(i)];
    }

    partition->y().setDataUpdated();
    partition->yp().setDataUpdated();

    // Evaluate this partition's residuals
    partition->evaluateResidual();

    auto* residual = partition->getResidual().getData();

    // Gather the residuals from the partition into the full monolithic vector
    for (size_t i = 0; i < partition->getInternalSize(); i++)
    {
      f[partition->getNodeConnection(i)] = residual[i];
    }
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
    ScaleMicrogridNetwork<ScalarT, IdxT>&          network,
    GridKit::PowerElectronicsModel<ScalarT, IdxT>& sys_model)
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