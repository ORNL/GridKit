#pragma once

#include <cassert>
#include <cmath>
#include <cstddef>
#include <functional>
#include <memory>
#include <optional>
#include <vector>

#include <sundials/sundials_types.h>

#include <GridKit/Model/PartitionEvaluator.hpp>
#include <GridKit/Solver/Dynamic/Interpolation.hpp>
#include <GridKit/Solver/Dynamic/Rosenbrock.hpp>

#include <resolve/Common.hpp>
#include <resolve/MemoryUtils.hpp>
#include <resolve/SystemSolver.hpp>
#include <resolve/matrix/Csr.hpp>
#include <resolve/matrix/MatrixHandler.hpp>
#include <resolve/vector/Vector.hpp>
#include <resolve/vector/VectorHandler.hpp>
#include <resolve/workspace/LinAlgWorkspace.hpp>
#include <resolve/workspace/LinAlgWorkspaceCpu.hpp>

namespace Integrator
{

  template <class ScalarT, typename IdxT>
  class MultiStageCosimulation
  {
    using State              = ReSolve::vector::Vector;
    using RealT              = typename GridKit::ScalarTraits<ScalarT>::RealT;
    using Rosenbrock         = Integrator::Rosenbrock<ScalarT, IdxT>;
    using PartitionEvaluator = GridKit::Model::PartitionEvaluator<ScalarT, IdxT>;
    using MemorySpace        = ReSolve::memory::MemorySpace;

  public:
    MultiStageCosimulation(std::vector<PartitionEvaluator*> partitions,
                           size_t                           num_stages     = 2,
                           size_t                           max_iterations = 2,
                           MemorySpace                      memspace       = MemorySpace::HOST);

    int initializeSimulation(RealT t0);

    int runSimulation(RealT tn, IdxT nout, std::optional<std::function<void(RealT)>> step_callback = {});

    int partitionStageSolve(int, double, double);

    int timestepper(const std::vector<double>&,
                    std::optional<std::function<void(double)>>);

    int distributeLocal(const State&);

    int distributeCoupling(const State&, int);

    int allocate();

  private:
    /**
     * @brief
     *
     */
    std::vector<std::unique_ptr<ReSolve::LinAlgWorkspaceCpu>> linear_workspaces_;

    /**
     * @brief initial time
     */
    RealT t0_;

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
     * @brief Number of stages used by the method
     */
    size_t num_states_;

    /**
     *
     * @brief The global state of the system
     */
    State y_;

    Rosenbrock::Parameters params_;

    Integrator::AdaptiveStep stepController_;

    /**
     *
     * @brief Memory space used
     */
    ReSolve::memory::MemorySpace memspace_;

    /**
     *
     * @brief Handles vector operations
     *
     */
    ReSolve::VectorHandler vector_handler_;

    /**
     *
     * @brief Handles matrix linear algebra operations
     *
     */
    ReSolve::MatrixHandler matrix_handler_;

    /**
     *
     * @brief linear solver for the sub-integrator
     *
     */
    std::vector<std::unique_ptr<ReSolve::SystemSolver>> lin_solver_;

    /**
     *
     * @brief Contains vectors of initial conditions for every step
     *
     */
    std::vector<std::unique_ptr<State>> x0_local_;

    /**
     *
     * @brief Vectors of coupling matrices for each partition
     *
     */
    std::vector<std::vector<std::unique_ptr<State>>> coupling_mat_;

    /**
     *
     * @brief partition solution multi-dimensional vectors
     *
     */
    std::vector<std::unique_ptr<State>> partition_solution_;

    /**
     *
     * @brief partition solution multi-dimensional vectors
     *
     */
    std::vector<ReSolve::matrix::Csr> partition_external_csr_;

    /**
     *
     * @brief partition solution multi-dimensional vectors
     *
     */
    std::vector<std::unique_ptr<ReSolve::matrix::Csr>> partition_internal_csr_;

    /**
     *
     *@brief Vector of partitions
     */
    std::vector<PartitionEvaluator*> partitions_;

    /**
     *
     *@brief Vector of solvers
     */
    std::vector<std::unique_ptr<Rosenbrock>> solvers_;
  };
} // namespace Integrator
