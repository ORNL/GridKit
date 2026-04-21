#pragma once

#include <cassert>
#include <cmath>
#include <cstddef>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include <sundials/sundials_types.h>

#include <GridKit/Model/PartitionEvaluator.hpp>
#include <GridKit/Solver/Dynamic/Rosenbrock.hpp>

#include <resolve/Common.hpp>
#include <resolve/MemoryUtils.hpp>
#include <resolve/SystemSolver.hpp>
#include <resolve/matrix/Csr.hpp>
#include <resolve/vector/Vector.hpp>
#include <resolve/vector/VectorHandler.hpp>
#include <resolve/workspace/LinAlgWorkspace.hpp>

namespace Integrator
{

  template <class ScalarT, typename IdxT>
  class MultiStageCosimulation
  {
    using State              = ReSolve::vector::Vector;
    using RealT              = typename GridKit::ScalarTraits<ScalarT>::RealT;
    using Rosenbrock         = Integrator::Rosenbrock<ScalarT, IdxT>;
    using PartitionEvaluator = GridKit::Model::PartitionEvaluator<ScalarT, IdxT>;

  public:
    int initializeSimulation(RealT t0);

    int runSimulation(RealT tn, IdxT nout, std::optional<std::function<void(RealT)>> step_callback = {});

    int partitionSolve();

    int timestepper(const std::vector<double>&                 out_times,
                    std::optional<std::function<void(double)>> out_cb = {});

    int distributeLocal(const State* global_y);

    int distributeCoupling(const State* global_y, int stage);

    int allocate();

  private:
    /**
     * @brief initial time
     */
    RealT  t0_;
    /**
     * @brief Maximum number of iterations per stage
     */
    size_t max_iterations_;

    /**
     *
     * @brief Number of partitions in the problem
     */
    size_t num_partitions_;

    /**
     *
     * @brief Number of stages used by the method
     */
    size_t num_stages_;

    /**
     *
     * @brief Memory space used
     */
    ReSolve::memory::MemorySpace memspace_;

    /**
     *
     * @brief Handles for linear algebra operations
     *
     */
    ReSolve::VectorHandler vector_handler_;

    /**
     *
     * @brief linear solver for the sub-integrator
     *
     */
    std::vector<ReSolve::SystemSolver> lin_solver;

    /**
     *
     * @brief Containes vectors of initial conditions for every step
     *
     */
    std::unique_ptr<std::unique_ptr<State>[]> x0_local_;

    /**
     *
     * @brief Vectors of coupling matrices for each partition
     *
     */
    std::unique_ptr<std::unique_ptr<State>[]> coupling_Mat_;

    /**
     *
     * @brief Vector of stages
     *
     */
    std::unique_ptr<State> stages_;

    /**
     *
     * @brief Vector of stages
     *
     */
    std::unique_ptr<State> partition_solution_;

    /**
     *
     *@brief Vector of partitions
     */
    std::vector<PartitionEvaluator*> partitions_;

    /**
     *
     *@brief Vector of solvers
     */
    std::vector<Rosenbrock*> solvers;
  }
} // namespace Integrator
