#pragma once

#include <cassert>
#include <cmath>
#include <cstddef>
#include <functional>
#include <memory>
#include <optional>
#include <vector>

#include <sundials/sundials_types.h>

#include <GridKit/LinearAlgebra/Solver/LinearSolver.hpp>
#include <GridKit/LinearAlgebra/Solver/ResolveSystemSolver.hpp>
#include <GridKit/LinearAlgebra/SparseMatrix/CsrMatrix.hpp>
#include <GridKit/LinearAlgebra/Vector/Vector.hpp>
#include <GridKit/LinearAlgebra/Vector/VectorHandler.hpp>
#include <GridKit/MemoryUtilities/MemoryUtils.hpp>
#include <GridKit/Model/Evaluator.hpp>
#include <GridKit/Model/PowerElectronics/SubsystemModel.hpp>
#include <GridKit/Model/PowerElectronics/SystemModelPowerElectronics.hpp>
#include <GridKit/ScalarTraits.hpp>
#include <GridKit/Solver/Dynamic/Ida.hpp>
#include <GridKit/Solver/Dynamic/Native/AdaptiveStep.hpp>
#include <GridKit/Solver/Dynamic/Native/ErrorNorm.hpp>
#include <GridKit/Solver/Dynamic/Native/InfNorm.hpp>
#include <GridKit/Solver/Dynamic/Native/StepController.hpp>
#include <GridKit/Solver/Dynamic/Rosenbrock.hpp>

#include "GridKit/Solver/Dynamic/Interpolation.hpp"
#include <resolve/workspace/LinAlgWorkspaceCpu.hpp>

namespace AnalysisManager
{

  template <class ScalarT, typename IdxT>
  class MultiStageCosimulation
  {

    using State          = GridKit::LinearAlgebra::Vector<ScalarT, IdxT>;
    using RealT          = typename GridKit::ScalarTraits<ScalarT>::RealT;
    using Interpolant    = GridKit::LagrangeInterpolant<ScalarT, IdxT>;
    using System         = GridKit::PowerElectronicsModel<ScalarT, IdxT>;
    using Subsystem      = GridKit::SubsystemModel<ScalarT, IdxT>;
    using MemorySpace    = GridKit::memory::MemorySpace;
    using CouplingMatrix = std::vector<std::shared_ptr<State>>;
    using IdaSolver      = AnalysisManager::Sundials::Ida<ScalarT, IdxT>;

  public:
    MultiStageCosimulation(
        size_t      max_iterations,
        size_t      num_stages,
        RealT       tolerance,
        MemorySpace memspace       = MemorySpace::HOST);

    ~MultiStageCosimulation();

    int allocate();

    int addSubsystems(Subsystem* subsys);

    int initializeSimulation(std::vector<ScalarT>& time_points, State* y, State* yp);

    int runSimulation(std::optional<std::function<void(const RealT t, const State&)>> step_callback = {});

    int zeroOrderHold(
        RealT t0,
        RealT dt);

    State& y()
    {
      return *y_;
    }

    State& yp()
    {
      return *yp_;
    }

  private:
    int partitionStageSolve(size_t curr_stage, size_t next_stage, RealT t0, RealT dt);

    int distributeLocalVariables(const State& global_y, const State& global_yp, Subsystem* partition);

    int gatherLocalVariables(State& global_y, State& global_yp, Subsystem* partition);

    int distributeCouplingVariables(size_t partition_idx, size_t stage);

    int distributeCouplingStage(size_t partition_idx, size_t stage_idx);

    int integrate(
        IdaSolver*                  solver,
        [[maybe_unused]] Subsystem* model,
        const std::vector<RealT>&   output_times,
        const bool                  is_last_iteration);

    /**
     * @brief Communication and output times for the simulation.
     */
    std::vector<RealT> time_points_;

    /**
     * @brief Number of stages used by the co-simulation method.
     */
    size_t num_stages_;

    /**
     * @brief Maximum number of coupling iterations performed at each stage.
     */
    size_t max_iterations_;

    /**
     * @brief Number of subsystem partitions of the problem.
     */
    size_t num_partitions_{0};

    /**
     * @brief Current global state of the complete system.
     */
    std::shared_ptr<State> y_;

    /**
     * @brief Current global state derivative of the complete system.
     *
     */
    std::shared_ptr<State> yp_;

    /**
     * @brief Global system states associated with current co-simulation stages.
     */
    std::vector<std::shared_ptr<State>> curr_stages_y_;

    /**
     * @brief Global system states associated with next co-simulation stages.
     *
     */
    std::vector<std::shared_ptr<State>> next_stages_y_;

    /**
     * @brief Global system states derivatives associated with the current co-simulation stages.
     *
     */
    std::vector<std::shared_ptr<State>> curr_stages_yp_;

    /**
     * @brief Global system states derivatives associated with next co-simulation stages.
     */
    std::vector<std::shared_ptr<State>> next_stages_yp_;

    /**
     * @brief Memory space used for  operations.
     */
    MemorySpace memspace_;

    /**
     * @brief Handler used to perform vector operations.
     */
    GridKit::LinearAlgebra::VectorHandler<ScalarT, IdxT> vector_handler_;

    /**
     * @brief Coupling state data associated with each subsystem partition.
     *
     */
    std::vector<CouplingMatrix> coupling_mat_y_;

    /**
     * @brief Coupling state derivative data associated with each subsystem partition.
     *
     */
    std::vector<CouplingMatrix> coupling_mat_yp_;

    /**
     * @brief Subsystem partitions that make up the complete system.
     */
    std::vector<Subsystem*> partitions_;

    /**
     * @brief IDA integrators associated with each subsystem partition.
     */
    std::vector<std::unique_ptr<IdaSolver>> ida_solvers_;

    /**
     * @brief Interpolants associated with each partition and stage on the on current states.
     */
    std::vector<std::vector<std::unique_ptr<Interpolant>>> interpolants_y_;

    /**
     * @brief Interpolants associated with each partition and stage on the on current states derivatives.
     */
    std::vector<std::vector<std::unique_ptr<Interpolant>>> interpolants_yp_;

    /**
     * @brief IDA sub-integrator tolerance.
     */
    RealT tol_;
  };
} // namespace AnalysisManager
