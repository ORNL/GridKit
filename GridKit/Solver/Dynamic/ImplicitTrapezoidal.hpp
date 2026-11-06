#pragma once

#include <cstddef>
#include <functional>
#include <memory>
#include <optional>
#include <vector>

#include <GridKit/LinearAlgebra/Solver/LinearSolver.hpp>
#include <GridKit/LinearAlgebra/Vector/Vector.hpp>
#include <GridKit/LinearAlgebra/Vector/VectorHandler.hpp>
#include <GridKit/MemoryUtilities/MemoryUtils.hpp>
#include <GridKit/Model/Evaluator.hpp>
#include <GridKit/ScalarTraits.hpp>
#include <GridKit/Solver/Dynamic/Native/ErrorNorm.hpp>

namespace AnalysisManager
{
  namespace NativeDynamicSolver
  {
    /**
     * @brief Fixed-step implicit trapezoidal integrator for index-one ODEs and DAEs.
     *
     * Nonlinear systems are solved with Newton's method and an Armijo backtracking
     * line search.
     */
    template <class ScalarT, typename IdxT>
    class ImplicitTrapezoidal
    {
      using RealT = typename GridKit::ScalarTraits<ScalarT>::RealT;
      using State = GridKit::LinearAlgebra::Vector<ScalarT, IdxT>;

    public:
      /**
       * @brief Parameters controlling fixed stepping and the nonlinear solve.
       */
      struct Parameters
      {
        /** @brief Fixed internal step size. */
        RealT  step_size_               = 1.0 / 240.0; // Quarter cycle
        /** @brief Maximum number of accepted internal steps. */
        size_t max_steps_               = 10000;
        /** @brief Maximum number of linear solves allowed in one nonlinear step. */
        size_t max_newton_iterations_   = 16;
        /** @brief Sufficient-decrease constant used by the Armijo line search. */
        RealT  armijo_constant_         = 1e-4;
        /** @brief Factor by which the line-search step length is reduced. */
        RealT  backtrack_factor_        = 0.5;
        /** @brief Smallest line-search step length that may be attempted. */
        RealT  minimum_step_length_     = 1e-8;
        /** @brief Residual-merit reduction ratio above which a Newton iteration is stagnant. */
        RealT  stagnation_ratio_        = 0.9;
        /** @brief Consecutive stagnant iterations that trigger a Jacobian refresh. */
        size_t max_stagnant_iterations_ = 2;
        /** @brief Maximum Jacobian refreshes permitted during one nonlinear step. */
        size_t max_jacobian_refreshes_  = 1;
      };

      /**
       * @brief Information about one accepted internal step.
       */
      struct StepInfo
      {
        /** @brief Simulation time at the end of the step. */
        RealT  sim_time_{};
        /** @brief Step size used for the step. */
        RealT  step_size_{};
        /** @brief Accepted step number, starting at one. */
        size_t step_no_{};
        /** @brief Number of Newton linear solves performed during the step. */
        size_t newton_iterations_{};
        /** @brief Number of rejected line-search trial points during the step. */
        size_t backtracks_{};
      };

      /**
       * @brief Running solver statistics since the last initializeSimulation() call.
       */
      struct Stats
      {
        /** @brief Number of accepted internal steps. */
        size_t num_steps_                = 0;
        /** @brief Number of Newton linear solves. */
        size_t num_newton_iterations_    = 0;
        /** @brief Number of model residual evaluations. */
        size_t num_residual_evaluations_ = 0;
        /** @brief Number of model Jacobian evaluations. */
        size_t num_jacobian_evaluations_ = 0;
        /** @brief Number of steps that reused a factorization from an earlier step. */
        size_t num_jacobian_reuses_      = 0;
        /** @brief Number of within-step Jacobian refreshes caused by poor convergence. */
        size_t num_jacobian_refreshes_   = 0;
        /** @brief Number of solves with a factored Newton matrix. */
        size_t num_linear_solves_        = 0;
        /** @brief Number of rejected line-search trial points. */
        size_t num_backtracks_           = 0;
        /** @brief Number of internal steps that failed to converge. */
        size_t num_convergence_failures_ = 0;
      };

      /**
       * @brief Construct an implicit trapezoidal integrator.
       *
       * The supplied model, linear solver, vector handler, and error norm must
       * remain valid for the lifetime of the integrator.
       *
       * @param[in,out] model Model defining \f$F(t,y,\dot{y})=0\f$.
       * @param[in,out] linear_solver Solver used for Newton systems.
       * @param[in,out] vector_handler Vector operations used by the integrator.
       * @param[in] error_norm Norm used to test Newton corrections for convergence.
       * @param[in] memspace Memory space in which vectors are stored and operated on.
       */
      ImplicitTrapezoidal(GridKit::Model::Evaluator<ScalarT, IdxT>*             model,
                          GridKit::LinearAlgebra::LinearSolver<ScalarT, IdxT>&  linear_solver,
                          GridKit::LinearAlgebra::VectorHandler<ScalarT, IdxT>& vector_handler,
                          const ErrorNorm<ScalarT, IdxT>*                       error_norm,
                          GridKit::memory::MemorySpace                          memspace = GridKit::memory::HOST);

      /**
       * @brief Allocate the integrator's internal vectors.
       *
       * @return Zero on success and a nonzero error code on failure.
       */
      [[nodiscard("May fail. Check error code.")]]
      int allocate();

      /**
       * @brief Initialize integration from the state currently exposed by the model.
       *
       * The model's `y()` and `yp()` vectors provide \f$y_0\f$ and
       * \f$\dot{y}_0\f$, respectively. This call also configures the linear
       * solver for the model's analytic CSR Jacobian and resets all statistics
       * and cached-factorization state.
       *
       * @param[in] t0 Initial simulation time.
       * @return Zero on success and a nonzero error code on failure.
       */
      [[nodiscard("May fail. Check error code.")]]
      int initializeSimulation(RealT t0);

      /**
       * @brief Integrate through a sequence of requested output times.
       *
       * Internal steps always use `Parameters::step_size_`. Output times between
       * accepted steps are reported using linear interpolation without changing
       * the accepted internal state. During an output callback, the model exposes
       * the state, derivative, and time associated with that output.
       *
       * @param[in] output_times Times at which output is requested, in increasing order.
       * @param[in] parameters Fixed-step and nonlinear-solver parameters.
       * @param[in] output_callback Optional callback invoked at every output time.
       * @param[in] step_callback Optional callback invoked after every accepted internal step.
       * @return Zero on success and a nonzero error code on failure.
       */
      [[nodiscard("May fail. Check error code.")]]
      int integrate(const std::vector<RealT>&                           output_times,
                    Parameters                                          parameters      = {},
                    std::optional<std::function<void(RealT)>>           output_callback = {},
                    std::optional<std::function<void(const StepInfo&)>> step_callback   = {});

      /**
       * @brief Advance the solution by one fixed step.
       *
       * Newton's method starts from the predictor
       * \f$y_{n+1}^{(0)}=y_n+dt\,\dot{y}_n\f$. If the nonlinear solve fails, the
       * last accepted model state and time are restored and the cached Jacobian
       * factorization is invalidated.
       *
       * @param[in] t0 Current accepted simulation time.
       * @param[in] dt Step size to take.
       * @param[in] parameters Nonlinear-solver parameters.
       * @return Zero on success and a nonzero error code on failure.
       */
      [[nodiscard("May fail. Check error code.")]]
      int timeStep(RealT t0, RealT dt, const Parameters& parameters = {});

      /** @brief Return running solver statistics. */
      const Stats& getStats() const
      {
        return stats_;
      }

      /** @brief Return the time of the current accepted internal state. */
      RealT getCurrentTime() const
      {
        return current_time_;
      }

    private:
      bool  validParameters(const Parameters& parameters) const;
      int   evaluateStepResidual(RealT t1, RealT dt, State& residual);
      int   evaluateStepJacobian(RealT t1, RealT dt);
      int   restoreAcceptedState();
      RealT residualMerit(State& residual) const;

      GridKit::Model::Evaluator<ScalarT, IdxT>*             model_;
      GridKit::LinearAlgebra::LinearSolver<ScalarT, IdxT>&  linear_solver_;
      GridKit::LinearAlgebra::VectorHandler<ScalarT, IdxT>& vector_handler_;
      const ErrorNorm<ScalarT, IdxT>*                       error_norm_;
      GridKit::memory::MemorySpace                          memspace_;

      RealT current_time_{};
      RealT previous_time_{};
      bool  initialized_{};
      bool  jacobian_factorized_{};
      bool  jacobian_valid_{};
      RealT jacobian_step_size_{};

      Stats    stats_;
      StepInfo last_step_info_;

      std::unique_ptr<State> current_state_;
      std::unique_ptr<State> current_derivative_;
      std::unique_ptr<State> previous_state_;
      std::unique_ptr<State> previous_derivative_;
      std::unique_ptr<State> iterate_;
      std::unique_ptr<State> residual_;
      std::unique_ptr<State> correction_;
    };
  } // namespace NativeDynamicSolver
} // namespace AnalysisManager
