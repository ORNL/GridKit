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

    using State                = GridKit::LinearAlgebra::Vector<ScalarT, IdxT>;
    using RealT                = typename GridKit::ScalarTraits<ScalarT>::RealT;
    using RosenbrockIntegrator = AnalysisManager::NativeDynamicSolver::Rosenbrock<ScalarT, IdxT>;
    using Interpolant          = GridKit::LagrangeInterpolant<ScalarT, IdxT>;
    using GridKitResolveSolver = GridKit::LinearAlgebra::ResolveSystemSolver<ScalarT, IdxT>;
    using Subsystem            = GridKit::SubsystemModel<ScalarT, IdxT>;
    using MemorySpace          = GridKit::memory::MemorySpace;
    using CouplingMatrix       = std::vector<std::shared_ptr<State>>;
    // using ErrorNorm            = NativeDynamicSolver::InfNorm<ScalarT, IdxT>;
    using IdaSolver            = AnalysisManager::Sundials::Ida<ScalarT, IdxT>;

    struct MergedTimePoints
    {
      std::vector<ScalarT> times;
      std::vector<bool>    from_first;
      std::vector<bool>    from_second;
    };

  public:
    MultiStageCosimulation(
        size_t      max_iterations = 2,
        size_t      num_stages     = 1,
        MemorySpace memspace       = MemorySpace::HOST);

    ~MultiStageCosimulation();

    int allocate();

    int addSubsystems(Subsystem* subsys);

    int initializeSimulation(std::vector<ScalarT>& time_points, State* y, State* yp);

    int runSimulation(std::optional<std::function<void(const State&)>> step_callback = {});

    State& y()
    {
      return y_;
    }

    void printStages() const
    {
      std::cout << "\nCurrent stages:\n";

      for (IdxT i = 0; i < y_.getSize(); ++i)
      {
        for (const auto& stage : curr_stages_y_)
        {
          const auto* data = stage->getData(memspace_);

          std::cout << std::setw(18)
                    << std::setprecision(10)
                    << std::scientific
                    << data[i];
        }

        std::cout << '\n';
      }

      std::cout << "\nNext stages:\n";

      for (IdxT i = 0; i < y_.getSize(); ++i)
      {
        for (const auto& stage : next_stages_y_)
        {
          const auto* data = stage->getData(memspace_);

          std::cout << std::setw(18)
                    << std::setprecision(10)
                    << std::scientific
                    << data[i];
        }

        std::cout << '\n';
      }

      std::cout << '\n';
    }

    void printState(std::string msg, State& y) const
    {
      std::cout << "\n"
                << msg << std::endl;

      for (IdxT i = 0; i < y.getSize(); ++i)
      {
        const auto* data = y.getData(memspace_);

        std::cout << std::setw(18)
                  << std::setprecision(10)
                  << std::scientific
                  << data[i];

        std::cout << '\n';
      }

      std::cout << '\n';
    }

    struct ResidualInfo
    {
      RealT  max_abs{0};
      size_t local_idx{0};
      IdxT   global_idx{0};
    };

    /**
     * @brief Evaluate F(t, y, y') for a partition with its current state and
     *        external data, and return the largest |F_i| and where it occurs.
     */
    ResidualInfo computeMaxResidual(Subsystem* partition, RealT t)
    {
      partition->updateTime(t, 0.0);
      partition->evaluateResidual();

      const auto* r = partition->getResidual().getData();

      ResidualInfo info;
      for (size_t i = 0; i < partition->getInternalSize(); ++i)
      {
        const RealT v = std::abs(r[i]);
        if (v > info.max_abs)
        {
          info.max_abs   = v;
          info.local_idx = i;
        }
      }
      info.global_idx = partition->getNodeConnection(info.local_idx);

      return info;
    }

    template <class ModelT>
    void printMaxResidual(const char* label, ModelT* model)
    {
      const auto* residual = model->getResidual().getData();

      double max_residual = 0.0;
      size_t max_index    = 0;

      for (size_t i = 0; i < model->getInternalSize(); ++i)
      {
        const double value = std::abs(residual[i]);

        if (value > max_residual)
        {
          max_residual = value;
          max_index    = i;
        }
      }

      std::cout << std::setprecision(16)
                << '\n'
                << label << '\n'
                << "  max residual = " << max_residual << '\n'
                << "  equation     = " << max_index << '\n';
    }

  private:
    int partitionStageSolve(size_t curr_stage, size_t next_stage, RealT t0, RealT dt);

    int distributeLocalVariables(const State& global_y, const State& global_yp, Subsystem* partition);

    int gatherLocalVariables(State& global_y, State& global_yp, Subsystem* partition);

    int distributeCouplingVariables(size_t p, size_t stage);

    int distributeCouplingStage(size_t partition_idx, size_t stage_idx);

    int ros_integrate(
        Subsystem*                        partition,
        RealT                             t0,
        RealT                             dt,
        const std::vector<RealT>&         output_times,
        const std::function<void(RealT)>& sol_out);

    int ZOH(RealT t0, RealT dt);

    MergedTimePoints mergeTimePoints(const std::vector<ScalarT>& first, const std::vector<ScalarT>& second);

    // int integrate(const std::vector<RealT>, std::function<void(RealT)> solout, Subsystem* model);

    int integrate(
        IdaSolver*                 solver,
        Subsystem*                 model,
        RealT                      t0,
        const std::vector<RealT>&  output_times,
        std::function<void(RealT)> sol_out);

    /**
     * @brief Linear algebra workspaces used by the partition solvers.
     */
    std::vector<std::unique_ptr<ReSolve::LinAlgWorkspaceCpu>> linear_workspaces_;

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
    State y_;

    /**
     * @brief Current global state derivative of the complete system.
     *
     */
    State yp_;

    /**
     * @brief Global system states associated with current co-simulation stages.
     */
    std::vector<std::unique_ptr<State>> curr_stages_y_;

    /**
     * @brief Global system states associated with next co-simulation stages.
     *
     */
    std::vector<std::unique_ptr<State>> next_stages_y_;

    /**
     * @brief Global system states associated with the co-simulation stages.
     *
     */
    std::vector<std::unique_ptr<State>> curr_stages_yp_;

    /**
     * @brief Global system states derivatives associated with next co-simulation stages.
     */
    std::vector<std::unique_ptr<State>> next_stages_yp_;

    /**
     * @brief Parameters used by the Rosenbrock partition integrators.
     */
    typename RosenbrockIntegrator::Parameters params_;

    /**
     * @brief Memory space used for vector storage and operations.
     */
    MemorySpace memspace_;

    /**
     * @brief Handler used to perform vector operations.
     */
    GridKit::LinearAlgebra::VectorHandler<ScalarT, IdxT> vector_handler_;

    /**
     * @brief Linear solver associated with each subsystem partition.
     *
     * Each solver uses the corresponding entry in linear_workspaces_.ReSolve::SystemSolver
     */
    std::vector<std::unique_ptr<ReSolve::SystemSolver>> resolve_lin_solver_;

    /** @brief Linear solver associated with each subsystem partition.
     *
     * Each solver uses the corresponding entry in linear_workspaces_.
     */
    std::vector<std::unique_ptr<GridKitResolveSolver>> gridkit_lin_solver_;

    /**
     * @brief Coupling data associated with each subsystem partition.
     *
     */
    std::vector<CouplingMatrix> coupling_mat_y_;

    /**
     * @brief Coupling data associated with each subsystem partition.
     *
     */
    std::vector<CouplingMatrix> coupling_mat_yp_;

    /**
     * @brief Subsystem partitions that make up the complete system.
     */
    std::vector<Subsystem*> partitions_;

    /**
     * @brief Rosenbrock integrators associated with each subsystem partition.
     *
     * Each integrator advances one partition using its corresponding
     * linear solver and coupling data.
     */
    std::vector<std::unique_ptr<RosenbrockIntegrator>> solvers_;

    /**
     * @brief IDA integrators associated with each subsystem partition.
     */
    std::vector<std::unique_ptr<IdaSolver>> ida_solvers_;
  };
} // namespace AnalysisManager
