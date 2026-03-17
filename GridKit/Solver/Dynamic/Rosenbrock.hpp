#pragma once

#include <cassert>
#include <cmath>
#include <cstddef>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include <GridKit/Model/Evaluator.hpp>

#include <resolve/Common.hpp>
#include <resolve/MemoryUtils.hpp>
#include <resolve/SystemSolver.hpp>
#include <resolve/matrix/Csr.hpp>
#include <resolve/vector/Vector.hpp>
#include <resolve/vector/VectorHandler.hpp>
#include <resolve/workspace/LinAlgWorkspace.hpp>

namespace Integrator
{
  using State = ReSolve::vector::Vector;

  /**
   * @brief Define control flow for `StepController`s to be able to control the step size of a `Rosenbrock` integrator.
   *
   */
  struct StepControl
  {
    /**
     * @brief Whether or not the step is accepted. A rejected step will cause the time step controller to discard
     *        the next state and re-step with the new `step_size`.
     *
     */
    bool accept;

    /**
     * @brief The step size the next step should take.
     *
     */
    double step_size;
  };

  /**
   * @brief Interface for step size controllers. Used by `Rosenbrock` integrators to decide when to accept/reject steps and
   *        what size each step should be.
   *
   */
  class StepController
  {
  public:
    /**
     * @brief Decide the control flow for the next step, based on information gathered by the integrator about the current step.
     *
     * @param err The estimated error made by the current step. Only calculated if `usesError()` returns true, otherwise it is assumed this
     *            value is not used.
     * @param prev_step The control flow from before the current step was taken. Can be used to accurately update the step size.
     * @param method_order The order of the method being used.
     * @return StepControl
     */
    virtual StepControl nextStep(double err, StepControl prev_step, uint8_t method_order) = 0;

    /**
     * @brief Return whether or not the `nextStep` method implementation uses the `err` parameter. If `false`, this parameter is not calculated.
     *
     */
    virtual bool usesError() const = 0;
  };

  /**
   * @brief Interface for error norms. Used to calculate the `err` parameter in `StepController::nextStep` based on a residual state error vector.
   *
   */
  class ErrorNorm
  {
  public:
    /**
     * @brief Calculate an error to be used by a step controller. Typically, an error > 1 indicates an error which does not meet tolerances, while
     *        an error < 1 indicates an error which meets tolerances. For that reason, tolerances should be included in the calculation of the error.
     *
     * @param err The state error residual being measured.
     * @param y The state that the error was calculated from. Can be used for proper relative error normalization.
     * @param yprev The state from the previous step. Can be used for proper relative error normalization.
     * @param handler A vector handler which can be used to facilitate vector operations.
     * @param memspace The memory space which vector operations should be performed in/.
     * @return double The error.
     */
    virtual double errorNorm(State& err, State& y, State& yprev, ReSolve::VectorHandler& handler, ReSolve::memory::MemorySpace memspace) const = 0;
  };

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
      double sim_time;

      /**
       * @brief The size of the step.
       *
       */
      double step_size;

      /**
       * @brief The size of the next step, as governed by the current `StepController` in use.
       *
       */
      double next_step_size;

      /**
       * @brief The estimated error made by the step, as calculated by the current `ErrorNorm` in use.
       *
       */
      double err_est;

      /**
       * @brief The step number, starting at 1.
       *
       */
      size_t step_no;

      /**
       * @brief Whether or not the integrator decided to skip computing the decomposition of the Jacobian on this step.
       *
       */
      bool skip_lu;

      /**
       * @brief Whether or not the integrator decided to skip evaluating the residual on the first stage on this step.
       *
       */
      bool skip_f;

      /**
       * @brief Whether or not this step was accepted by the `StepController` in use.
       *
       */
      bool accepted;

      std::string csv_report() const;
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
      std::vector<StepInfo> rejections;

      /**
       * @brief Information of each step which the integrator decided to skip re-factoring the Jacobian.
       *
       */
      std::vector<StepInfo> skip_lu_steps;

      /**
       * @brief How many steps the integrator has taken.
       *
       */
      size_t num_steps = 0;

      /**
       * @brief Number of model residual function evaluations.
       *
       */
      size_t f_evals = 0;

      /**
       * @brief Number of model residual function evaluations which have been skipped by the integrator.
       *
       */
      size_t f_skipped = 0;

      /**
       * @brief Number of model Jacobian evaluations.
       *
       */
      size_t jac_evals = 0;

      /**
       * @brief Number of linear solves against the model Jacobian.
       *
       */
      size_t decomp_solves = 0;

      /**
       * @brief Minimum step size.
       *
       */
      double min_step = INFINITY;

      /**
       * @brief Maximum step size.
       *
       */
      double max_step = 0;

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
      double starting_step = 1e-5;

      /**
       * @brief The maximum number of steps the integrator should take. If the integrator has not reached the final time before
       *        taking this many steps, then integration is stopped. For more details, see `integrate()`.
       *
       */
      size_t max_steps = 2000;

      /**
       * @brief Whether or not the integrator should attempt to skip Jacobian decompositions.
       *
       * @note This feature is only available if the underlying method is a Rosenbrock-W method and can speed up
       *       the time taken to compute each step. However, the overall number of steps taken will increase.
       *
       */
      bool skip_lu = false;
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
      size_t num_stages;

      /**
       * @brief The coefficient along the diagonal of the Gamma matrix.
       *
       */
      RealT gamma;

      /**
       * @brief A vector of sums of rows of the alpha matrix. These are the classic
       *        Runge-Kutta 'c' coefficients, or abscissae. The size of this vector
       *        should be equal to `num_stages`.
       *
       */
      std::unique_ptr<RealT[]> alpha_sum;

      /**
       * @brief A vector of sums of rows of the Gamma matrix. The size of this vector
       *        should be equal to `num_stages`.
       *
       */
      std::unique_ptr<RealT[]> gamma_sum;

      /**
       * @brief A vector of weights for constructing the final solution from the stages.
       *        The size of this vector should be equal to `num_stages`.
       *
       */
      std::unique_ptr<RealT[]> m;

      /**
       * @brief OPTIONAL vector of coefficients for the embedded error method. If it exists,
       *        the size of this vector should be equal to `num_stages`.
       *
       */
      std::unique_ptr<RealT[]> e;

      /**
       * @brief The transformed A coefficient matrix. Strictly lower triangular and stored in dense row-major form.
       *        Upper triangular terms are not accessed. Should be `num_stages` by `num_stages` large.
       *
       */
      std::unique_ptr<RealT[]> A;

      /**
       * @brief The transformed C coefficient matrix. Strictly lower triangular and stored in dense row-major form.
       *        Upper triangular terms are not accessed. Should be `num_stages` by `num_stages` large.
       *
       */
      std::unique_ptr<RealT[]> C;

      /**
       * @brief OPTIONAL matrix of dense coefficients. Defines how the stages should be transformed into interpolant
       *        nodes for computing dense output. The interpolating polynomial has an order one less than the order of
       *        the method, and two interpolant nodes are already pre-computed, so if this matrix exists it should be
       *        `order` - 2 by `num_stages` large.
       *
       */
      std::unique_ptr<RealT[]> H;

      /**
       * @brief What ODE order these coefficients satisfy. If `is_dae` is true, then the coefficients must additionally satisfy
       *        DAE conditions up to this order. If `is_w` is true, then the coefficients must additionally satisfy ROW conditions
       *        up to this order. If `is_krylov` is true, then the coefficients must additionally satisfy ROK condition up to this order.
       *
       */

      uint8_t order;

      /**
       * @brief Whether or not these coefficients are appropriate to use in a Rosenbrock-Krylov (ROK) solver.
       *
       */
      bool is_krylov;

      /**
       * @brief Whether or not these coefficients satisfy Rosenbrock-W (ROW) order conditions up to `order`.
       *        The integrator may take advantage of this fact by e.g. using time-delay Jacobians to speed up computation.
       *
       */
      bool is_w;

      /**
       * @brief Whether or not these coefficients satisfy DAE order conditions up to `order`. If this is not true,
       *        these coefficients should not be used to solve models with algebraic conditions (indicated by a
       *        `Model::Evaluator::tag_` value of 0).
       *
       */
      bool is_dae;

      /**
       * @brief Whether or not this tableau contains an embedded error estimator method.
       *
       */
      constexpr bool has_embedded() const
      {
        return static_cast<bool>(e);
      }

      /**
       * @brief Whether or not this tableau contains coefficients which can be used to construct dense output.
       *
       */
      constexpr bool hasDenseOutput() const
      {
        return static_cast<bool>(H);
      }

      /**
       * @brief Helper function for accessing elements of `A`
       *
       */
      constexpr RealT getA(size_t row, size_t col) const
      {
        return A[row * num_stages + col];
      }

      constexpr bool                     can_reuse_asum(size_t stage) const;
      constexpr bool                     can_reuse_asum_for_out() const;
      constexpr std::tuple<bool, size_t> error_estimator_stage() const;

      static Tableau lin_implicit_euler();
      static Tableau rodas5p();
    };

  private:
    /**
     * @brief The current step size of the integrator. If integration is ever stopped and resumed, this value will be used for
     *        the initial step after resuming.
     *
     */
    double step_size_ = 0;

    /**
     * @brief The step size of the previous step. Used for operations which need to be done on the current step, but step size
     *        control for the next step has already been performed.
     *
     */
    double prev_step_size_ = 0;

    /**
     * @brief Whether or not the integrator should attempt to skip Jacobian decomposition on the next step. Controlled by the
     *        time stepping algorithm in `integrate()`. Generally, this should only be set if we suspect the Jacobian for the
     *        last step is a good enough approximation of the Jacobian on the next step. Non-ROW methods need exact Jacobians,
     *        so this should only be set for ROW methods, and when the step size for the next step is the same as the step size
     *        as the previous step.
     *
     */
    bool skip_lu_ = false;

    /**
     * @brief Whether or not the integrator should attempt to skip the residual function evaluation of the first stage on the
     *        next step. This should only be used when a step is rejected and the residual function is evaluated at the exact
     *        same arguments as the previous step. Then `RHS_first_stage_` can be re-used rather than re-calculated.
     *
     */
    bool skip_f_ = false;

    /**
     * @brief Keeps track of whether or not the integrator currently has valid dense coefficients.
     *        i.e. they have been computed and haven't been invalidated by taking another step. This can be used to avoid
     *        re-computing dense coefficients when interpolating states multiple times in one step.
     *
     */
    bool dense_coefficients_valid_ = false;

    /**
     * @brief The tableau of Rosenbrock coefficients currently being used by the integrator.
     *
     */
    Tableau tab_;

    /**
     * @brief The model being simulated.
     *
     */
    GridKit::Model::Evaluator<ScalarT, IdxT>* model_;

    /**
     * @brief The linear solver to be used during integration in `time_step()`.
     *
     */
    ReSolve::SystemSolver&  lin_solver_;
    /**
     * @brief The vector handler to be used for vector operations by the integrator.
     *
     */
    ReSolve::VectorHandler& vector_handler_;

    /**
     * @brief The `ErrorNorm` to be used by the `StepController` in `integrate()`.
     *
     * @note Can be `nullptr` if no `ErrorNorm` is configured, in case the `StepController` does not need an error to be calculated.
     *
     */
    const ErrorNorm* err_norm_;

    /**
     * @brief The memory space where linear algebra operations hsould be done in.
     *
     */
    ReSolve::memory::MemorySpace memspace_;

    /**
     * @brief The current simulation time.
     *
     */
    double current_time_;

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
     * @brief Running statistics of the integrator. Reset by `initializeSimulation()`.
     *
     */
    Stats stats_;

  public:
    Rosenbrock(Tableau&&                                 tab,
               GridKit::Model::Evaluator<ScalarT, IdxT>* model,
               ReSolve::SystemSolver&                    lin_solver,
               ReSolve::VectorHandler&                   vector_handler,
               const ErrorNorm*                          err_norm,
               ReSolve::memory::MemorySpace              memspace = ReSolve::memory::HOST);

    [[nodiscard("May fail. Check error code.")]]
    int allocate();

    [[nodiscard("May fail. Check error code.")]]
    int initializeSimulation(RealT t0);

    [[nodiscard("May fail. Check error code.")]]
    int integrate(const std::vector<double>&                          out_times,
                  StepController&                                     step_controller,
                  Parameters                                          params  = {},
                  std::optional<std::function<void(double)>>          out_cb  = {},
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
       *        Stored separately from `RHS_` since it can be re-used by future steps if the current step is rejected.
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
       * @brief Jacobian matrix
       *
       */
      std::unique_ptr<ReSolve::matrix::Csr> jacobian_;

      /**
       * @brief Whether or not the Jacobian has been factorized. Initial factorization is required,
       *        but pivots can be re-used for a refactorization later. Re-factorization can introduce
       *        errors when the Jacobian differs significantly, so this flag should be reset when
       *        a significant event can happen, such as re-initializing the simulation.
       *
       * @see `initializeSimulation()`
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
     * @see `calc_dense_coeff()`, `interp_dense()`
     *
     */
    std::unique_ptr<std::unique_ptr<State>[]> dense_coeff_;

  public:
    [[nodiscard("May fail. Check error code.")]]
    int time_step(double t0, double dt);

    State& error_estimate() const;

    [[nodiscard("May fail. Check error code.")]]
    int calc_dense_coeff();

    [[nodiscard("May fail. Check error code.")]]
    int interp_dense(double theta);
  };

  /**
   * @brief A simple textbook adaptive `StepController` which seeks to meet a relative and absolute tolerance
   *        based on an error estimate.
   *
   */
  class AdaptiveStep : public StepController
  {
    /**
     * @brief Parameters for the step controller.
     *
     */
    struct Parameters
    {
      /**
       * @brief The minimum multiple by which the step size can be multiplied to obtain the new step size.
       *        Increasing this can allow the integrator to be slightly more conservative in selecting the step size
       *        - decreasing the number of steps taken but increasing the risk of failing the next step.
       *
       * @note Should be between 0 and 1.
       *
       */
      double facmin = 0.2;

      /**
       * @brief The maximum multiple by which the step size can be multiplied to obtain the new step size.
       *        Decreasing this will make the integrator more conservative in selecting the step size -
       *        increasing the number of steps taken but decreasing the risk of failing the next step.
       *
       * @note Should be greater than 1.
       *
       */
      double facmax = 5.0;

      /**
       * @brief A "fudge factor" introduced to decrease risk of failing a step. The larger the fudge factor,
       *        the more likely steps will fail, but fewer steps will be taken.
       *
       * @note Should be between 0 and 1.
       *
       */
      double facscale = 0.9;
    } params_;

  public:
    AdaptiveStep(const Parameters& params) : params_(params)
    {
    }

    StepControl nextStep(double err, StepControl prev_step, uint8_t method_order) final;

    /**
     * @brief This controller uses error estimates.
     *
     * @see `nextStep()`
     *
     */
    constexpr bool usesError() const final
    {
      return true;
    }
  };

  /**
   * @brief A fixed step controller which doesn't change the step size and accepts every step.
   *        Useful if you know what time scale your simulation operates on apriori and you're
   *        using a method without an embedded error controller.
   *
   *        To set the fixed size, set the `Rosenbrock::Parameters::starting_step` parameter.
   *
   */
  class FixedStep : public StepController
  {
    StepControl nextStep(double err, StepControl prev_step, uint8_t method_order) final;

    /**
     * @brief This controller does not use error estimates.
     *
     * @see `nextStep()`
     *
     */
    constexpr bool usesError() const final
    {
      return false;
    }
  };

  class InfNorm : public ErrorNorm
  {
    mutable struct
    {
      std::unique_ptr<State> out_;
      std::unique_ptr<State> scale_;
      std::unique_ptr<State> yprev_abs_;
    } workspace_;

  public:
    struct Parameters
    {
      std::unique_ptr<State> atol;
      double                 rtol;
    } params_;

    InfNorm(Parameters&& params) : params_(std::move(params))
    {
    }

    double errorNorm(State& err, State& y, State& yprev, ReSolve::VectorHandler& handler, ReSolve::memory::MemorySpace memspace) const final;
  };
} // namespace Integrator
