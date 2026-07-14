#pragma once

#include <cassert>
#include <cmath>
#include <cstddef>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include <GridKit/LinearAlgebra/Solver/LinearSolver.hpp>
#include <GridKit/LinearAlgebra/Vector/Vector.hpp>
#include <GridKit/LinearAlgebra/Vector/VectorHandler.hpp>
#include <GridKit/MemoryUtilities/MemoryUtils.hpp>
#include <GridKit/Model/Evaluator.hpp>
#include <GridKit/ScalarTraits.hpp>
#include <GridKit/Solver/Dynamic/Native/ErrorNorm.hpp>
#include <GridKit/Solver/Dynamic/Native/StepController.hpp>

namespace Integrator
{
  /**
   * @brief A Rosenbrock integrator designed to simulate a `GridKit::Model::Evaluator` model.
   *        Rosenbrock integrators are best for quick medium-accuracy simulation and
   *        boast a large amount of customization to simulating a given model.
   *
   *        For the list of available Rosenbrock methods, see `Rosenbrock::Tableau`.
   */
  template <class ScalarT, typename IdxT>
  class Rosenbrock
  {
    using RealT = typename GridKit::ScalarTraits<ScalarT>::RealT;
    using State = GridKit::LinearAlgebra::Vector<ScalarT, IdxT>;

  public:
    /**
     * @brief Keeps track of a variety of notable properties of a single step.
     *
     */
    struct StepInfo
    {
      /**
       * @brief The simulation time at the beginning of the step
       *
       */
      RealT  sim_time_;
      /**
       * @brief The size of the step.
       *
       */
      RealT  step_size_;
      /**
       * @brief The size of the next step, as governed by the current `StepController` in use.
       *
       */
      RealT  next_step_size_;
      /**
       * @brief The estimated error made by the step, as calculated by the current `ErrorNorm` in use.
       *
       */
      RealT  err_est_;
      /**
       * @brief The step number, starting at 1.
       *
       */
      size_t step_no_;
      /**
       * @brief Whether or not the integrator decided to skip computing the decomposition of the Jacobian on this step.
       *
       */
      bool   skip_lu_;
      /**
       * @brief Whether or not the integrator decided to skip evaluating the residual on the first stage on this step.
       *
       */
      bool   skip_f_;
      /**
       * @brief Whether or not this step was accepted by the `StepController` in use.
       *
       */
      bool   accepted_;

      std::string csvReport() const;
      std::string report() const;
    };

    /**
     * @brief Keeps track of running statistics of the integrator since the last `initializeSimulation()` call.
     *
     */
    struct Stats
    {
      /**
       * @brief Information of each step which has been rejected.
       *
       */
      std::vector<StepInfo> rejections_;
      /**
       * @brief Information of each step which the integrator decided to skip re-factoring the Jacobian.
       *
       */
      std::vector<StepInfo> skip_lu_steps_;
      /**
       * @brief How many steps the integrator has taken.
       *
       */
      size_t                num_steps_     = 0;
      /**
       * @brief Number of model residual function evaluations.
       *
       */
      size_t                f_evals_       = 0;
      /**
       * @brief Number of model residual function evaluations which have been skipped by the integrator.
       *
       */
      size_t                f_skipped_     = 0;
      /**
       * @brief Number of model Jacobian evaluations.
       *
       */
      size_t                jac_evals_     = 0;
      /**
       * @brief Number of linear solves against the model Jacobian.
       *
       */
      size_t                decomp_solves_ = 0;
      /**
       * @brief Minimum step size.
       *
       */
      RealT                 min_step_      = INFINITY;
      /**
       * @brief Maximum step size.
       *
       */
      RealT                 max_step_      = 0;

      std::string report() const;
      Stats&      operator+=(const Stats& other);
    };

    /**
     * @brief Parameters of the integrator.
     *
     */
    struct Parameters
    {
      /**
       * @brief What step size the first step should take.
       *
       * @todo Consider adding a starting step size selector to select this automatically.
       */
      RealT  starting_step_ = 1e-5;
      /**
       * @brief The maximum number of steps the integrator should take. If the integrator has not reached the final time before
       *        taking this many steps, then integration is stopped. For more details, see `integrate()`.
       *
       */
      size_t max_steps_     = 2000;
      /**
       * @brief Whether or not the integrator should attempt to skip Jacobian decompositions.
       *
       * @note This feature is only available if the underlying method is a Rosenbrock-W method and can speed up
       *       the time taken to compute each step. However, the overall number of steps taken will increase.
       *
       */
      bool   skip_lu_       = false;
    };

    /**
     * @brief A list of the coefficients needed to complete Rosenbrock integration in an accurate way.
     *        Each tableau can be considered a different method or scheme, and will have different properties,
     *        advantages, and disadvantages.
     *
     */
    struct Tableau
    {
      /**
       * @brief The number of stages used by the method. Each stage requires one model residual evaluation.
       *
       */
      size_t                   num_stages_;
      /**
       * @brief The coefficient along the diagonal of the Gamma matrix.
       *
       */
      RealT                    gamma_;
      /**
       * @brief A vector of sums of rows of the alpha matrix. These are the classic
       *        Runge-Kutta 'c' coefficients, or abscissae. The size of this vector
       *        should be equal to `num_stages`.
       *
       */
      std::unique_ptr<RealT[]> alpha_sum_;
      /**
       * @brief A vector of sums of rows of the Gamma matrix. The size of this vector
       *        should be equal to `num_stages`.
       *
       */
      std::unique_ptr<RealT[]> gamma_sum_;
      /**
       * @brief A vector of weights for constructing the final solution from the stages.
       *        The size of this vector should be equal to `num_stages`.
       *
       */
      std::unique_ptr<RealT[]> m_;
      /**
       * @brief OPTIONAL vector of coefficients for the embedded error method. If it exists,
       *        the size of this vector should be equal to `num_stages`.
       *
       */
      std::unique_ptr<RealT[]> e_;
      /**
       * @brief The transformed A coefficient matrix. Strictly lower triangular and stored in dense row-major form.
       *        Upper triangular terms are not accessed. Should be `num_stages` by `num_stages` large.
       *
       */
      std::unique_ptr<RealT[]> A_;
      /**
       * @brief The transformed C coefficient matrix. Strictly lower triangular and stored in dense row-major form.
       *        Upper triangular terms are not accessed. Should be `num_stages` by `num_stages` large.
       *
       */
      std::unique_ptr<RealT[]> C_;
      /**
       * @brief OPTIONAL matrix of dense coefficients. Defines how the stages should be transformed into interpolant
       *        nodes for computing dense output. The interpolating polynomial has an order one less than the order of
       *        the method, and two interpolant nodes are already pre-computed, so if this matrix exists it should be
       *        `order` - 2 by `num_stages` large.
       *
       */
      std::unique_ptr<RealT[]> H_;
      /**
       * @brief What ODE order these coefficients satisfy. If `is_dae` is true, then the coefficients must additionally satisfy
       *        DAE conditions up to this order. If `is_w` is true, then the coefficients must additionally satisfy ROW conditions
       *        up to this order. If `is_krylov` is true, then the coefficients must additionally satisfy ROK condition up to this order.
       *
       */

      uint8_t order_;
      /**
       * @brief Whether or not these coefficients are appropriate to use in a Rosenbrock-Krylov (ROK) solver.
       *
       */
      bool    is_krylov_;
      /**
       * @brief Whether or not these coefficients satisfy Rosenbrock-W (ROW) order conditions up to `order`.
       *        The integrator may take advantage of this fact by e.g. using time-delay Jacobians to speed up computation.
       *
       */
      bool    is_w_;
      /**
       * @brief Whether or not these coefficients satisfy DAE order conditions up to `order`. If this is not true,
       *        these coefficients should not be used to solve models with algebraic conditions (indicated by a
       *        `Model::Evaluator::tag_` value of 0).
       *
       */
      bool    is_dae_;

      /**
       * @brief Whether or not this tableau contains an embedded error estimator method.
       *
       */
      constexpr bool hasEmbedded() const
      {
        return static_cast<bool>(e_);
      }

      /**
       * @brief Whether or not this tableau contains coefficients which can be used to construct dense output.
       *
       */
      constexpr bool hasDenseOutput() const
      {
        return static_cast<bool>(H_);
      }

      /**
       * @brief Helper function for accessing elements of `A`
       *
       */
      constexpr RealT getA(size_t row, size_t col) const
      {
        return A_[row * num_stages_ + col];
      }

      constexpr bool                  canReuseAsum(size_t stage) const;
      constexpr bool                  canReuseAsumForOut() const;
      constexpr std::optional<size_t> errorEstimatorStage() const;

      static Tableau linImplicitEuler();
      static Tableau rodas5p();
    };

  private:
    /**
     * @brief The current step size of the integrator. If integration is ever stopped and resumed, this value will be used for
     *        the initial step after resuming.
     *
     */
    RealT step_size_                = 0;
    /**
     * @brief The step size of the previous step. Used for operations which need to be done on the current step, but step size
     *        control for the next step has already been performed.
     *
     */
    RealT prev_step_size_           = 0;
    /**
     * @brief Whether or not the integrator should attempt to skip Jacobian decomposition on the next step. Controlled by the
     *        time stepping algorithm in \link integrate() \endlink . Generally, this should only be set if we suspect the Jacobian for the
     *        last step is a good enough approximation of the Jacobian on the next step. Non-ROW methods need exact Jacobians,
     *        so this should only be set for ROW methods, and when the step size for the next step is the same as the step size
     *        as the previous step.
     *
     */
    bool  skip_lu_                  = false;
    /**
     * @brief Whether or not the integrator should attempt to skip the residual function evaluation of the first stage on the
     *        next step. This should only be used when a step is rejected and the residual function is evaluated at the exact
     *        same arguments as the previous step. Then \ref RHS_first_stage_ can be re-used rather than re-calculated.
     *
     */
    bool  skip_f_                   = false;
    /**
     * @brief Keeps track of whether or not the integrator currently has valid dense coefficients.
     *        i.e. they have been computed and haven't been invalidated by taking another step. This can be used to avoid
     *        re-computing dense coefficients when interpolating states multiple times in one step.
     *
     */
    bool  dense_coefficients_valid_ = false;

    /**
     * @brief The tableau of Rosenbrock coefficients currently being used by the integrator.
     *
     */
    Tableau                                               tab_;
    /**
     * @brief The model being simulated.
     *
     */
    GridKit::Model::Evaluator<ScalarT, IdxT>*             model_;
    /**
     * @brief The linear solver to be used during integration in \link timeStep() \endlink.
     *
     */
    GridKit::LinearAlgebra::LinearSolver<ScalarT, IdxT>&  lin_solver_;
    /**
     * @brief The vector handler to be used for vector operations by the integrator.
     *
     */
    GridKit::LinearAlgebra::VectorHandler<ScalarT, IdxT>& vector_handler_;
    /**
     * @brief The `ErrorNorm` to be used by the `StepController` in \link integrate() \endlink.
     *
     * @note Can be `nullptr` if no `ErrorNorm` is configured, in case the `StepController` does not need an error to be calculated.
     *
     * @todo Should be removed from the `Rosenbrock` class. Whether or not this is needed is dependent on the `StepController`, so it should be stored there.
     */
    const ErrorNorm<ScalarT, IdxT>*                       err_norm_;
    /**
     * @brief The memory space where linear algebra operations hsould be done in.
     *
     */
    GridKit::memory::MemorySpace                          memspace_;

    /**
     * @brief The current simulation time.
     *
     */
    RealT current_time_;

    /**
     * @brief The state from last step.
     *
     */
    std::unique_ptr<State> y_prev_;
    /**
     * @brief The state on the current step.
     *
     */
    std::unique_ptr<State> y_cur_;
    /**
     * @brief The incoming state for the next step.
     *
     */
    std::unique_ptr<State> y_new_;
    /**
     * @brief Used as output for dense output interpolation.
     *
     */
    std::unique_ptr<State> y_interp_;

    /**
     * @brief Configured parameters for the integrator.
     *
     */
    Parameters params_;

    /**
     * @brief Running statistics of the integrator. Reset by \link initializeSimulation() \endlink.
     *
     */
    Stats stats_;

  public:
    Rosenbrock(Tableau&&                                             tab,
               GridKit::Model::Evaluator<ScalarT, IdxT>*             model,
               GridKit::LinearAlgebra::LinearSolver<ScalarT, IdxT>&  lin_solver,
               GridKit::LinearAlgebra::VectorHandler<ScalarT, IdxT>& vector_handler,
               const ErrorNorm<ScalarT, IdxT>*                       err_norm,
               GridKit::memory::MemorySpace                          memspace = GridKit::memory::HOST);

    [[nodiscard("May fail. Check error code.")]]
    int allocate();

    [[nodiscard("May fail. Check error code.")]]
    int initializeSimulation(RealT t0);

    [[nodiscard("May fail. Check error code.")]]
    int integrate(const std::vector<RealT>&                           out_times,
                  StepController<RealT>&                              step_controller,
                  Parameters                                          params  = {},
                  std::optional<std::function<void(RealT)>>           out_cb  = {},
                  std::optional<std::function<void(const StepInfo&)>> step_cb = {});

  private:
    /**
     * @brief A workspace used by `time_step` for the necessary linear algebra computations.
     *
     */
    mutable struct
    {
      /**
       * @brief Sum of A-weighted states used for function evaluation in a stage.
       *
       */
      std::unique_ptr<State> asum_;
      /**
       * @brief Sum of C-weighted states used in linearization in a stage.
       *
       */
      std::unique_ptr<State> csum_;
      /**
       * @brief Right-hand side of linear solve used in a stage.
       *
       */
      std::unique_ptr<State> RHS_;
      /**
       * @brief Right-hand side of linear solve used in the first stage.
       *        Stored separately from \ref RHS_ since it can be re-used by future steps if the current step is rejected.
       *
       * @see `skip_f_`
       *
       */
      std::unique_ptr<State> RHS_first_stage_;
      /**
       * @brief Time Jacobian for non-autonomous models. Currently unused, since models are assumed to be autonomous.
       *
       */
      std::unique_ptr<State> dFdt_;
      /**
       * @brief Vector representation of a diagonal mass matrix.
       *
       */
      std::unique_ptr<State> mass_;
      /**
       * @brief Estimated error produced by a step in a method with an empbedded error estimator.
       *
       * @see `Tableau::e`
       *
       */
      std::unique_ptr<State> err_est_;

      /**
       * @brief Whether or not the Jacobian has been factorized. Initial factorization is required,
       *        but pivots can be re-used for a refactorization later. Re-factorization can introduce
       *        errors when the Jacobian differs significantly, so this flag should be reset when
       *        a significant event can happen, such as re-initializing the simulation.
       *
       * @see \link initializeSimulation() \endlink
       *
       */
      bool jacobian_analyzed_ = false;

      /**
       * @brief Integration stages - each a low-order approximation of the state at the abscissae given by
       *        the tableau.
       *
       */
      std::unique_ptr<std::unique_ptr<State>[]> stages_;
    } workspace_;

    /**
     * @brief Dense output interpolation nodes. Used to generate output states in-between steps.
     *
     * @see \link calcDenseCoeff() \endlink, \link interpDense() \endlink
     *
     */
    std::unique_ptr<std::unique_ptr<State>[]> dense_coeff_;

  public:
    [[nodiscard("May fail. Check error code.")]]
    int timeStep(RealT t0, RealT dt);

    State& errorEstimate() const;

    [[nodiscard("May fail. Check error code.")]]
    int calcDenseCoeff();

    [[nodiscard("May fail. Check error code.")]]
    int interpDense(RealT theta);
  };
} // namespace Integrator
