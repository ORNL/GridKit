// PartitionUtilities.hpp
#pragma once

#include <cassert>
#include <cstddef>
#include <unordered_map>
#include <vector>

#include <GridKit/Model/Evaluator.hpp>
#include <GridKit/Model/PowerElectronics/PartitionInterface/BusLinePartitionInterface.hpp>
#include <GridKit/Model/PowerElectronics/SubsystemModel.hpp>
#include <GridKit/ScalarTraits.hpp>
#include <GridKit/Testing/Testing.hpp>
#include <GridKit/Utilities/MapFromCsr.hpp>

#include <examples/PowerElectronics/ExamplesHelper/MicrogridNetwork.hpp>

/**
 * @brief Partition a scaled microgrid network  and assemble into subsystem models.
 *
 * Divides the network into contiguous groups of IBRs and creates one
 * subsystem for each group. The IBRs are distributed as evenly as possible,
 * with any remainder assigned to the leading partitions. A partition
 * interface is added the appropriate subsystem
 *
 * @param[in,out] network Microgrid network to be partitioned.
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
 *
 * @note The subsystem models are dynamically allocated. The caller is
 *       responsible for releasing and deleting them.
 */
template <class ScalarT, typename IdxT>
void partitionMicrogridNetwork(
    ScaleMicrogridNetwork<ScalarT, IdxT>&                                   network,
    std::vector<GridKit::PowerElectronics::SubsystemModel<ScalarT, IdxT>*>& subsystems,
    size_t                                                                  num_partitions)
{
  const size_t num_ibrs = 2 * network.N_size;

  auto interface_id = network.model_id_next;

  assert(num_partitions <= num_ibrs);

  subsystems.resize(num_partitions);

  IdxT q = num_ibrs / num_partitions;
  IdxT r = num_ibrs % num_partitions;

  IdxT partition_begin = 0;

  for (IdxT j = 0; j < num_partitions; j++)
  {
    auto* partition = new GridKit::PowerElectronics::SubsystemModel<ScalarT, IdxT>();

    // Add the reference signal node to the first partition.
    if (j == 0)
    {
      partition->addNode(&network.dg_signal);
    }

    // Divide the IBRs as evenly as possible among the partitions.
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

    // Insert a partition interface inbetween every partition
    if (partition_end < num_ibrs)
    {

      const IdxT boundary_line_index = partition_end - 1;

      auto* bus_interface = new GridKit::PowerElectronics::BusLinePartitionInterface<ScalarT, IdxT>(
          &network.buses[boundary_line_index],
          network.lines[boundary_line_index],
          interface_id++);

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
 *
 * @post Each subsystem residual has been evaluated at the supplied
 *       @p time and @p alpha.
 * @post @p f contains the reconstructed residual for all internal variables
 *       represented by the subsystems.
 */
template <class ScalarT, typename IdxT>
void parallelFunctionEval(
    const std::vector<GridKit::PowerElectronics::SubsystemModel<ScalarT, IdxT>*>& subsystems,
    const std::vector<ScalarT>&                                                   y,
    const std::vector<ScalarT>&                                                   yp,
    std::vector<ScalarT>&                                                         f,
    typename GridKit::ScalarTraits<ScalarT>::RealT                                time,
    typename GridKit::ScalarTraits<ScalarT>::RealT                                alpha)
{
#ifdef _OPENMP
#pragma omp parallel for schedule(guided)
#endif
  for (size_t p = 0; p < subsystems.size(); ++p)
  {
    auto* partition = subsystems[p];
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
 * @brief Verify that the subsystem Jacobian is an exact subset of the
 * full-system Jacobian to ensure accuracy of subsystem Jacobians.
 *
 * Each subsystem Jacobian entry should match the corresponding entry in the
 * full-system Jacobian after converting local subsystem indices back to global
 * indices. The check verifies that:
 *
 * - every subsystem entry exists in the full-system Jacobian,
 * - every full-system entry whose row and column belong to the subsystem
 * exists in the subsystem Jacobian, and
 * - matching entries agree within the specified tolerance.
 *
 * @note When the components in a subsystem are evaluated in a different order
 *       than in the monolithic reference system, a larger tolerance may be
 *       required to account for floating-point roundoff introduced by the
 *       different summation order.
 *
 * @param full_jac Full-system Jacobian.
 * @param subsystem Subsystem associated with sub_jac.
 * @param tolerance Maximum allowed difference between matching entries.
 *
 * @return true if the subsystem Jacobian is an exact subset of the
 * full-system Jacobian, false otherwise.
 */

template <typename ScalarT, typename IdxT>
bool verifySubsystemJacobian(
    GridKit::LinearAlgebra::CsrMatrix<ScalarT, IdxT>*         full_jac,
    GridKit::PowerElectronics::SubsystemModel<ScalarT, IdxT>& subsystem,
    ScalarT                                                   tolerance = std::numeric_limits<ScalarT>::epsilon())
{
  const auto subsystem_jacobian = GridKit::Testing::MapFromCsr(subsystem.getCsrJacobian());
  const auto reference_jacobian = GridKit::Testing::MapFromCsr(full_jac);

  const auto partition_size = subsystem.getInternalSize();

  // A table to lookup if global connection index falls within this partition
  std::unordered_map<IdxT, IdxT> global_to_local;

  for (size_t i = 0; i < partition_size; ++i)
  {
    global_to_local[subsystem.getNodeConnection(i)] = i;
  }

  for (size_t row = 0; row < partition_size; ++row)
  {
    const auto global_row = subsystem.getNodeConnection(row);

    const auto& reference_row = reference_jacobian[global_row];
    const auto& partition_row = subsystem_jacobian[row];

    size_t matched = 0;

    // Compare only full-system entries belonging to the subsystem.
    for (const auto& [global_col, reference_value] : reference_row)
    {
      // Skip columns outside the subsystem.
      if (!global_to_local.contains(global_col))
      {
        continue;
      }

      const auto local_col = global_to_local.at(global_col);

      // The corresponding subsystem entry must exist.
      if (!partition_row.contains(local_col))
      {
        return false;
      }

      // Compare the Jacobian values.
      bool status = GridKit::Testing::isEqual(partition_row.at(local_col), reference_value, tolerance);

      if (!status)
      {
        return false;
      }
      ++matched;
    }

    // Make sure there are no additional entries in the subsystem row that hasn't been checked.
    if (matched != partition_row.size())
    {
      return false;
    }
  }

  return true;
}
