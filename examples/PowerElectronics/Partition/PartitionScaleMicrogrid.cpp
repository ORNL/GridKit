#include <omp.h>

#include <cassert>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <string>

#include <GridKit/Model/PowerElectronics/PartitionInterface/BusPartitionInterface.hpp>
#include <GridKit/Model/PowerElectronics/SubsystemModel.hpp>
#include <GridKit/Solver/Dynamic/DynamicSolver.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>
#include <GridKit/Testing/Testing.hpp>

#include "PowerElectronicsExamplesHelper/JacTestHelper.hpp"
#include "PowerElectronicsExamplesHelper/MicrogridNetwork.hpp"
#include "PowerElectronicsExamplesHelper/PartitionUtilities.hpp"

using index_type = size_t;
using real_type  = double;

using SignalNode  = GridKit::PowerElectronics::SignalNode<real_type, index_type>;
using Bus         = GridKit::PowerElectronics::MicrogridBus<real_type, index_type>;
using BusDQ       = GridKit::MicrogridBusDQ<real_type, index_type>;
using DGGenerator = GridKit::DistributedGenerator<real_type, index_type>;
using Line        = GridKit::MicrogridLine<real_type, index_type>;
using Load        = GridKit::MicrogridLoad<real_type, index_type>;

/*
 * Output data for each parallel function evaluation run
 */
struct RunResult
{
  bool        success = true;
  index_type  num_partitions;
  real_type   partition_eval_time;  // seconds
  real_type   monolithic_eval_time; // seconds
  real_type   speedup;              // monolithic / partition
  real_type   max_error;
  std::string jacobian_status;
};

/*
 * Reference data from the monolithic system evaluation.
 *
 * The monolithic system is evaluated once and its states, residual, and
 * evaluation time are stored here. Each partition configuration reuses this
 * data for validation and performance comparison.
 */
struct MonolithicReference
{
  std::vector<real_type> y;
  std::vector<real_type> yp;
  std::vector<real_type> residual;
  real_type              monolithic_eval_time;
};

/**
 * @brief Evaluate the monolithic system and store reference data.
 *
 * Initializes the global state and state-derivative vectors with deterministic
 * values, evaluates and times the monolithic residual, evaluates the full
 * Jacobian, and stores the resulting data for later comparison with
 * partitioned evaluations.
 *
 * The residual evaluation time is measured independently from the Jacobian
 * evaluation so that the stored timing represents only the cost of the
 * monolithic residual evaluation.
 *
 * @param[in,out] system Allocated monolithic power electronics system.
 *
 * @return MonolithicReference containing the state vector, state-derivative
 *         vector, residual, and residual evaluation time.
 *
 * @pre @p system is non-null.
 * @pre @p system has already been assembled and allocated.
 *
 * @post The state and state-derivative vectors of @p system contain the
 *       deterministic reference values.
 * @post The monolithic residual and Jacobian have been evaluated.
 * @post The returned reference contains an independent copy of the
 *       monolithic residual.
 */
MonolithicReference evaluateMonolithicSystem(GridKit::PowerElectronicsModel<real_type, index_type>* system)
{
  MonolithicReference reference;

  // ---------------------------------------------------------------------------
  // Initialize the reference state and state-derivative vectors
  // ---------------------------------------------------------------------------
  reference.y.resize(system->size());
  reference.yp.resize(system->size());

  for (index_type i = 0; i < system->size(); ++i)
  {
    reference.y[i]  = static_cast<real_type>(i + 1);
    reference.yp[i] = static_cast<real_type>(i + 1);
  }

  // Copy the reference values into the vectors owned by the system model.
  auto* system_y  = system->y().getData();
  auto* system_yp = system->yp().getData();

  for (index_type i = 0; i < system->size(); ++i)
  {
    system_y[i]  = reference.y[i];
    system_yp[i] = reference.yp[i];
  }

  system->y().setDataUpdated();
  system->yp().setDataUpdated();

  // ---------------------------------------------------------------------------
  // Evaluate and time the monolithic residual
  // ---------------------------------------------------------------------------

  auto start_time = std::chrono::high_resolution_clock::now();

  system->updateTime(0.1, 0.1);
  system->evaluateResidual();

  auto end_time = std::chrono::high_resolution_clock::now();

  reference.monolithic_eval_time = std::chrono::duration<real_type>(end_time - start_time).count();

  // evaluate the monolithic Jacobian
  system->evaluateJacobian();

  // Store a copy of the residual return reference data
  auto* residual = system->getResidual().getData();

  reference.residual.assign(residual, residual + system->size());

  return reference;
}

/**
 * @brief Evaluate one partition configuration and validate it against the
 *        monolithic reference system.
 *
 * Creates the requested subsystem decomposition, allocates each subsystem,
 * evaluates the partitioned residual, verifies each subsystem Jacobian against
 * the monolithic Jacobian, and compares the reconstructed global residual with
 * the stored monolithic residual.
 *
 * The partitioned residual evaluation is timed independently and compared with
 * the previously measured monolithic residual evaluation time to compute the
 * resulting speedup.
 *
 * @param[in,out] network Scaled microgrid network to partition.
 * @param[in] system Allocated monolithic system used as the reference.
 * @param[in] reference Stored monolithic state, derivative, residual, and
 *                      timing information.
 * @param[in] num_partitions Number of subsystem partitions to create.
 *
 * @return Performance and validation results for the requested partition count.
 *
 * @pre @p system is non-null and has already been evaluated.
 * @pre @p network has been constructed using buildScaleMicrogridNetwork().
 * @pre @p num_partitions is greater than zero and does not exceed the number
 *      of IBRs in the network.
 * @pre @p reference contains state and derivative vectors consistent with the
 *      size of @p system.
 *
 * @post All temporary subsystem models created by this function are released
 *       and deleted before returning.
 * @post The returned result reports the partition timing, speedup, residual
 *       error, and subsystem Jacobian validation status.
 */
template <class ScalarT, typename IdxT>
RunResult evaluatePartitioning(
    ScaleMicrogridNetwork<ScalarT, IdxT>&                  network,
    GridKit::PowerElectronicsModel<real_type, index_type>* system,
    const MonolithicReference&                             reference,
    index_type                                             num_partitions)
{

  // ---------------------------------------------------------------------------
  // Create and allocate subsystem partitions
  // ---------------------------------------------------------------------------

  std::vector<GridKit::SubsystemModel<real_type, index_type>*> subsystems;

  partitionNetwork(network, subsystems, num_partitions);

  for (auto* partition : subsystems)
  {
    partition->allocate();
  }

  // Global residual reconstructed from the subsystem residuals.
  std::vector<real_type> f(system->size(), 1.0);
  // Elementwise error between the monolithic and reconstructed residuals.
  std::vector<real_type> error(system->size(), 1.0);

  // ---------------------------------------------------------------------------
  // Evaluate and time the partitioned residual
  // ---------------------------------------------------------------------------
  auto start_time = std::chrono::high_resolution_clock::now();

  evaluatePartitionResiduals(subsystems, reference.y, reference.yp, f, 0.1, 0.1);

  auto end_time = std::chrono::high_resolution_clock::now();

  auto partition_eval_time = std::chrono::duration<real_type>(end_time - start_time);

  // ---------------------------------------------------------------------------
  // Verify subsystem Jacobians against the monolithic Jacobian
  // ---------------------------------------------------------------------------

  auto* system_jacobian = system->getCsrJacobian();

  bool      jacobian_match = true;
  real_type Jac_tol        = 1e-13;

  for (auto* partition : subsystems)
  {
    partition->evaluateJacobian();

    jacobian_match = jacobian_match
                     && GridKit::Testing::verifySubsystemJacobian(*system_jacobian,
                                                                  *partition->getCsrJacobian(),
                                                                  *partition,
                                                                  Jac_tol);
  }

  // ---------------------------------------------------------------------------
  // Compare the reconstructed and monolithic residuals
  // ---------------------------------------------------------------------------
  real_type max_error = 0.0;

  for (index_type i = 0; i < system->size(); ++i)
  {
    error[i] = std::abs(f[i] - reference.residual[i]) / (std::abs(reference.residual[i]) + 1.0);

    if (max_error < error[i])
    {
      max_error = error[i];
    }
  }

  // ---------------------------------------------------------------------------
  // Store performance and validation results
  // ---------------------------------------------------------------------------
  RunResult result;

  result.num_partitions       = num_partitions;
  result.partition_eval_time  = partition_eval_time.count();
  result.monolithic_eval_time = reference.monolithic_eval_time;
  result.speedup              = reference.monolithic_eval_time / partition_eval_time.count();
  result.max_error            = max_error;
  result.jacobian_status      = jacobian_match ? "Correct" : "Wrong";

  // ---------------------------------------------------------------------------
  // Check validation results
  // ---------------------------------------------------------------------------
  if (!jacobian_match)
  {
    std::cout << "ERROR: At least one subsystem Jacobian is incorrect!"
              << std::endl;

    result.success = false;
  }

  if (max_error > std::numeric_limits<real_type>::epsilon())
  {
    std::cout << "ERROR: Max Error too high!: " << max_error << std::endl;
    result.success = false;
  }

  // ---------------------------------------------------------------------------
  // Release and destroy temporary subsystem models
  // ---------------------------------------------------------------------------

  for (auto* partition : subsystems)
  {
    partition->release();
    delete partition;
  }

  return result;
}

/**
 * @brief Print the header for the partition evaluation results table.
 */
void printResultsHeader()
{
  std::cout << std::left
            << std::setw(16) << "num_partitions"
            << std::right
            << std::setw(16) << "partition_time"
            << std::setw(18) << "monolithic_time"
            << std::setw(12) << "speedup"
            << std::setw(14) << "error"
            << std::setw(16) << "Jacobians"
            << '\n';

  std::cout << std::string(92, '-') << '\n';
}

/**
 * @brief Print one row of partition evaluation results.
 *
 * @param[in] result Results from the partitioned system evaluation.
 */
void printResult(const RunResult& result)
{
  std::cout << std::left
            << std::setw(16) << result.num_partitions
            << std::right
            << std::setw(14) << std::fixed << std::setprecision(4)
            << result.partition_eval_time << " s"
            << std::setw(16) << result.monolithic_eval_time << " s"
            << std::setw(11) << std::setprecision(2)
            << result.speedup << "x"
            << std::setw(14) << std::scientific << std::setprecision(3)
            << result.max_error << ' '
            << std::setw(16) << result.jacobian_status
            << '\n';

  std::cout << std::defaultfloat;
}

/**
 * @brief Benchmark and validate partitioned residual evaluation of a large
 *        scaled microgrid. It also confirms the correctness of the subsystem
 *        Jacobian.
 *
 * Builds a scaled microgrid once, constructs a monolithic reference system,
 * and compares several subsystem decompositions of the same physical network.
 *
 * For each requested partition count, the benchmark measures the partitioned
 * residual evaluation time, computes its speedup relative to the monolithic
 * evaluation, verifies the reconstructed residual, and checks every subsystem
 * Jacobian against the monolithic Jacobian.
 *
 * The physical network and monolithic system are constructed only once.
 * Individual subsystem models are created and destroyed separately for each
 * partition count.
 *
 * @return 0 if all residual and Jacobian validation checks pass; otherwise 1.
 */
int main(int argc, char const* argv[])
{
  index_type N_size = 5000;

  std::vector<index_type> num_partitions_list = {500};

  /*
   * If command-line arguments are provided, the first argument specifies
   * N_size and all remaining arguments specify partition counts.
   *
   * Example:
   *   ./PartitionScaleMicrogrid 5000 10 48 100 500
   */
  if (argc > 1)
  {
    try
    {
      N_size = static_cast<index_type>(std::stoull(argv[1]));

      if (N_size < 1)
      {
        std::cerr << "ERROR: N_size must be at least 1.\n";
        return 1;
      }

      // When N_size is supplied explicitly, at least one partition count
      // must also be supplied.
      if (argc < 3)
      {
        std::cerr << "ERROR: At least one partition count must be provided "
                  << "when N_size is specified.\n";
        return 1;
      }

      num_partitions_list.clear();

      for (int i = 2; i < argc; ++i)
      {
        index_type num_partitions = static_cast<index_type>(std::stoull(argv[i]));

        if (num_partitions < 1)
        {
          std::cerr << "ERROR: Number of partitions must be at least 1.\n";
          return 1;
        }

        if (num_partitions > 2 * N_size)
        {
          std::cerr << "ERROR: Number of partitions ("
                    << num_partitions
                    << ") cannot exceed the number of IBRs ("
                    << 2 * N_size
                    << ").\n";
          return 1;
        }

        num_partitions_list.push_back(num_partitions);
      }
    }
    catch (const std::exception& e)
    {
      std::cerr
          << "ERROR: Invalid command-line argument: "
          << e.what()
          << "\n";

      return 1;
    }
  }

  const bool use_jac = true;

  // Build the physical network once.
  ScaleMicrogridNetwork<real_type, index_type> network(N_size);

  // Build, assemble and allocate the monolithic system once.
  auto* system = new GridKit::PowerElectronicsModel<real_type, index_type>(use_jac);

  assembleSystemLeftToRight(network, *system);
  system->allocate();

  // Evaluate the monolithic reference once.
  const MonolithicReference reference = evaluateMonolithicSystem(system);

  printResultsHeader();

  // Only the partitioned system is rebuilt and evaluated for each partition count.
  for (const index_type p : num_partitions_list)
  {
    assert(p <= 2 * N_size);

    // Partition the network and perform the parallel function evaluation.
    const RunResult r = evaluatePartitioning(network, system, reference, p);

    if (!r.success)
    {
      delete system;
      return 1;
    }

    // Output the partition evaluation results.
    printResult(r);
  }

  delete system;

  return 0;
}