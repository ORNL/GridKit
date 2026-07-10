#include "Rosenbrock.hpp"

#include <format>
#include <iomanip>
#include <sstream>

#include <GridKit/Constants.hpp>
#include <GridKit/MemoryUtilities/ResolveMemoryUtils.hpp>

/**
 * @brief A small helper macro to "bubble" errors. The Rosenbrock implementations call many
 * fallible model and linear algebra functions, and often the only thing that can be
 * done when one of these fails is to fail the current function as well. This macro
 * wraps a fallible function call and causes the current function to fail as well,
 * returning the same error code as the called function.
 *
 */
#define BUBBLE_FAIL(arg) \
  do                     \
  {                      \
    if (int err = (arg)) \
      return err;        \
  } while (false)

namespace Integrator
{
  /**
   * @brief Constructs a single-line report string out of the `StepInfo` object which is suitable to be included in a CSV file.
   *
   * Useful if you would like to plot step info.
   *
   */
  template <class ScalarT, typename IdxT>
  std::string Rosenbrock<ScalarT, IdxT>::StepInfo::csvReport() const
  {
    std::stringstream out;
    out << std::scientific << std::setprecision(20)
        << sim_time_ << ','
        << step_size_ << ','
        << next_step_size_ << ','
        << err_est_ << ','
        << step_no_ << ','
        << skip_lu_ << ','
        << skip_f_ << ','
        << accepted_;

    return out.str();
  }

  /**
   * @brief Constructs a human-readable report string out of the `StepInfo` object.
   *
   * Useful for debug dumps.
   *
   */
  template <class ScalarT, typename IdxT>
  std::string Rosenbrock<ScalarT, IdxT>::StepInfo::report() const
  {
    std::stringstream out;
    out << std::scientific << std::setprecision(20)
        << "Simulation Time: " << sim_time_ << '\n'
        << "Step Size:       " << step_size_ << '\n'
        << "Next Step Size:  " << next_step_size_ << '\n'
        << "Error Estimate:  " << err_est_ << '\n'
        << "Step Number:     " << step_no_ << '\n'
        << "Skip LU:         " << skip_lu_ << '\n'
        << "Skip F:          " << skip_f_ << '\n'
        << "Accepted:        " << accepted_;

    return out.str();
  }

  /**
   * @brief Constructs a human-readable report string out of the `Stats` object.
   *
   * Useful for reporting at the end of a simulation.
   *
   */
  template <class ScalarT, typename IdxT>
  std::string Rosenbrock<ScalarT, IdxT>::Stats::report() const
  {
    std::stringstream out;
    out << "Rejections: " << rejections_.size()
        << "\nSteps: " << num_steps_
        << "\nSkip LU Steps: " << skip_lu_steps_.size()
        << "\nMin Step: " << min_step_
        << "\nMax Step: " << max_step_
        << "\nRHS Function Evals: " << f_evals_
        << "\nRHS Function Skipped: " << f_skipped_
        << "\nJacobian Evals: " << jac_evals_
        << "\nLinear Solves: " << decomp_solves_;

    return out.str();
  }

  /**
   * @brief Accumulate statistics into one object. Useful to collect statistics over runs of several simulations, such as when
   *        there is a discrete event.
   *
   * @param other The other statistics to add to this one.
   *
   * @todo Right now, the step numbers for \ref rejections_ and \ref skip_lu_steps_ are impossible to tell apart from the different simulations.
   */
  template <class ScalarT, typename IdxT>
  typename Rosenbrock<ScalarT, IdxT>::Stats& Rosenbrock<ScalarT, IdxT>::Stats::operator+=(const Stats& other)
  {
    rejections_.insert(rejections_.end(), other.rejections_.begin(), other.rejections_.end());
    skip_lu_steps_.insert(skip_lu_steps_.end(), other.skip_lu_steps_.begin(), other.skip_lu_steps_.end());

    num_steps_     += other.num_steps_;
    f_evals_       += other.f_evals_;
    f_skipped_     += other.f_skipped_;
    jac_evals_     += other.jac_evals_;
    decomp_solves_ += other.decomp_solves_;

    min_step_ = std::min(min_step_, other.min_step_);
    max_step_ = std::max(max_step_, other.max_step_);

    return *this;
  }

  /**
   * @brief Checks to see if the \ref asum_ variable can be re-used from the previous stage.
   *
   * Often, a row in the Rosenbrock A matrix is the exact same as the previous row, but with an additional
   * non-zero element at the end. In this case, the \ref asum_ variable is the exact same as the previous stage,
   * but with one extra additional term. The integrator can take advantage of this and reduce a matmul for computing
   * \ref asum_ down to a single `axpy`.
   *
   * Typically, \ref asum_ is only initialized if it needs to be re-calculated from scratch. For this reason, this
   * function will always return `false` for `stage == 0`, forcing \ref asum_ to be initialized for the first stage.
   *
   * @param stage The stage being checked.
   */
  template <class ScalarT, typename IdxT>
  constexpr bool Rosenbrock<ScalarT, IdxT>::Tableau::canReuseAsum(size_t stage) const
  {
    assert(stage < num_stages_);

    if (stage == 0)
      return false;
    else
    {
      for (size_t j = 0; j < stage - 1; j++)
      {
        if (getA(stage, j) != getA(stage - 1, j))
        {
          return false;
        }
      }
      return true;
    }
  }

  /**
   * @brief Checks to see if the \ref asum_ variable can be re-used from the last stage to compute the output state
   *
   * Often, the Rosenbrock m vector is the exact same as the last row in the A matrix, but with an additional
   * non-zero element at the end. In this case, the output state for a step is equal to the \ref asum_ variables from
   * the last stage but with a single weighted vector added to it. The integrator can take advantage of this and reduce
   * a matmul for computing the final state down to a single `axpy`.
   *
   * If a method has only a single stage, then \ref asum_ will not have been initialized (as it will just be equal to y0).
   * Therefore, this method will return `false`.
   *
   */
  template <class ScalarT, typename IdxT>
  constexpr bool Rosenbrock<ScalarT, IdxT>::Tableau::canReuseAsumForOut() const
  {
    if (num_stages_ == 1)
      return false;

    for (size_t j = 0; j < num_stages_ - 1; j++)
    {
      if (getA(num_stages_ - 1, j) != m_[j])
      {
        return false;
      }
    }

    return true;
  }

  /**
   * @brief Returns the index of a stage which can be used as an embedded error estimator, if that stage exists.
   *
   * Sometimes, Rosenbrock methods are designed in such a way where a stage can be used as an embedded error estimator.
   * Typically, the embedded error estimator is a linear combination of the stages, so in this particular case the
   * weights will be 1 on the estimator stage and 0 on all other stages.
   *
   * @pre This function is only valid to call if there is an embedded error estimator and its coefficients
   * are included in this tableau.
   */
  template <class ScalarT, typename IdxT>
  constexpr std::optional<size_t> Rosenbrock<ScalarT, IdxT>::Tableau::errorEstimatorStage() const
  {
    assert(e_);

    std::optional<size_t> re;
    for (size_t j = 0; j < num_stages_; j++)
    {
      if (e_[j] == GridKit::ONE<RealT> && !re)
      {
        re = j;
      }
      else if (e_[j] != GridKit::ZERO<RealT>)
      {
        return {};
      }
    }
    return re;
  }

  /**
   * @brief Construct a new Rosenbrock integrator.
   *
   * @param tab The tableau to be used for this integrator. Since tableaus contain `std::unique_ptr`, it must be moved into the integrator.
   * @param model The model to be simulated. Despite taking a pointer, this must be a valid pointer to an `Evaluator`.
   * Must have \ref GridKit::Model::Evaluator::tag() set, and must be in Hessenberg form (\f(F(\dot y, y) = \dot y - f(y)\f)).
   * @param lin_solver The linear solver to be used when constructing stages during simulation. The reference must remain valid for as long
   * as the Rosenbrock integrator lives.
   * @param vector_handler The vector handler to be used when simulating. The reference must remain valid for as long as the Rosenbrock
   * integrator lives.
   * @param err_norm The error norm to use in calculating error for the `StepController`. Will not be accessed if the `StepController`
   * does not need error (such as `FixedStep`), so `nullptr` can be passed in that circumstance.
   * @param memspace The memory space that linear algebra operations should be performed in.
   */
  template <class ScalarT, typename IdxT>
  Rosenbrock<ScalarT, IdxT>::Rosenbrock(Tableau&&                                             tab,
                                        GridKit::Model::Evaluator<ScalarT, IdxT>*             model,
                                        ReSolve::SystemSolver&                                lin_solver,
                                        GridKit::LinearAlgebra::VectorHandler<ScalarT, IdxT>& vector_handler,
                                        const ErrorNorm<ScalarT, IdxT>*                       err_norm,
                                        GridKit::memory::MemorySpace                          memspace)
    : tab_(std::move(tab)),
      model_(model),
      lin_solver_(lin_solver),
      vector_handler_(vector_handler),
      err_norm_(err_norm),
      memspace_(memspace)
  {
  }

  /**
   * @brief Allocates memory based on the the size of \ref model_. Must be called before other methods.
   *
   * @note This method can fail.
   *
   * @pre Allocates the Jacobian matrix, which requires knowledge of the number of nonzero elements. This must be known at this point,
   * so `allocate()` must be called beforehand on \ref model_ to count the nonzero elements and allocate the CSR Jacobian.
   *
   * @return An error code, with 0 as success.
   */
  template <class ScalarT, typename IdxT>
  int Rosenbrock<ScalarT, IdxT>::allocate()
  {
    size_t size = static_cast<size_t>(model_->size());

    y_prev_                     = std::make_unique<State>(size);
    y_cur_                      = std::make_unique<State>(size);
    y_new_                      = std::make_unique<State>(size);
    y_interp_                   = std::make_unique<State>(size);
    workspace_.asum_            = std::make_unique<State>(size);
    workspace_.csum_            = std::make_unique<State>(size);
    workspace_.RHS_             = std::make_unique<State>(size);
    workspace_.RHS_first_stage_ = std::make_unique<State>(size);
    workspace_.dFdt_            = std::make_unique<State>(size);
    workspace_.mass_            = std::make_unique<State>(size);

    resolve_rhs_ = std::make_unique<ReSolve::vector::Vector>(size);
    resolve_lhs_ = std::make_unique<ReSolve::vector::Vector>(size);

    BUBBLE_FAIL(y_prev_->allocate(memspace_));
    BUBBLE_FAIL(y_cur_->allocate(memspace_));
    BUBBLE_FAIL(y_new_->allocate(memspace_));
    BUBBLE_FAIL(y_interp_->allocate(memspace_));
    BUBBLE_FAIL(workspace_.asum_->allocate(memspace_));
    BUBBLE_FAIL(workspace_.csum_->allocate(memspace_));
    BUBBLE_FAIL(workspace_.RHS_->allocate(memspace_));
    BUBBLE_FAIL(workspace_.RHS_first_stage_->allocate(memspace_));
    BUBBLE_FAIL(workspace_.dFdt_->allocate(memspace_));
    BUBBLE_FAIL(workspace_.mass_->allocate(memspace_));

    if (tab_.e_)
    {
      std::optional<size_t> err_est_stage = tab_.errorEstimatorStage();
      if (!err_est_stage)
      {
        workspace_.err_est_ = std::make_unique<State>(size);
        BUBBLE_FAIL(workspace_.err_est_->allocate(memspace_));
      }
    }

    workspace_.stages_ = std::make_unique<std::unique_ptr<State>[]>(tab_.num_stages_);
    for (size_t i = 0; i < tab_.num_stages_; i++)
    {
      workspace_.stages_[i] = std::make_unique<State>(size);
      BUBBLE_FAIL(workspace_.stages_[i]->allocate(memspace_));
      workspace_.stages_[i]->setToZero(memspace_);
    }

    if (tab_.order_ > 2)
    {
      dense_coeff_ = std::make_unique<std::unique_ptr<State>[]>(tab_.order_ - 2);
      for (size_t i = 0; i < static_cast<size_t>(tab_.order_ - 2); i++)
      {
        dense_coeff_[i] = std::make_unique<State>(size);
        BUBBLE_FAIL(dense_coeff_[i]->allocate(memspace_));
      }
    }

    workspace_.jacobian_ = std::make_unique<ReSolve::matrix::Csr>(size, size, model_->getCsrJacobian()->getNnz());

    return 0;
  }

  /**
   * @brief Initializes the simulation. Must be called before \ref integrate() or \ref timeStep(). Must also be called after
   * discrete events.
   *
   * - Sets the simulation time to `t0` and copies the initial condition from \ref model_.
   * - Analyzes \ref model_ Jacobian sparsity and runs the preconditioner
   * - Generates the mass matrix from \ref GridKit::Model::Evaluator::tag(). If the tag is not properly set, then initialization will fail.
   * - Resets \ref stats_.
   *
   * @note This method can fail.
   *
   * @pre Must have called \ref allocate(). The `model.tag_` variable must be properly constructed.
   *
   * @todo Document mass matrix construction
   *
   * @param t0 The starting simulation time.
   * @return An error code, with 0 as success.
   */
  template <class ScalarT, typename IdxT>
  int Rosenbrock<ScalarT, IdxT>::initializeSimulation(RealT t0)
  {
    current_time_ = t0;
    BUBBLE_FAIL(y_cur_->copyFromExternal(model_->y().data(), memspace_, memspace_));
    workspace_.jacobian_analyzed_ = false;

    GridKit::LinearAlgebra::CsrMatrix<RealT, IdxT>* model_jacobian = model_->getCsrJacobian();
    BUBBLE_FAIL(workspace_.jacobian_->setDataPointers(
        model_jacobian->getRowData(),
        model_jacobian->getColData(),
        model_jacobian->getValues(),
        GridKit::memory::memorySpaceAsResolve(memspace_)));
    BUBBLE_FAIL(lin_solver_.setMatrix(workspace_.jacobian_.get()));
    BUBBLE_FAIL(lin_solver_.analyze());

    // TODO: This function needs to be called to properly use a preconditioner in ReSolve (if there is any), but currently will error
    // if there is no preconditioner configured. Once we can detect if a preconditioner is configured, we can restore this functionality.
    // Also, we will always want to use *right* preconditioning.
    // BUBBLE_FAIL(lin_solver_.preconditionerSetup("right"));

    if (model_->tag().size() != static_cast<size_t>(model_->size()))
    {
      std::cerr << "Model tag is either unset or does not match the size of the model\n";
      return 1;
    }

    std::unique_ptr<RealT[]> mass = std::make_unique<RealT[]>(model_->tag().size());
    for (size_t i = 0; i < static_cast<size_t>(model_->size()); i++)
    {
      mass[i] = model_->tag()[i] ? GridKit::ONE<RealT> : GridKit::ZERO<RealT>;
    }
    BUBBLE_FAIL(workspace_.mass_->copyFromExternal(mass.get(), memspace_, memspace_));

    stats_ = Stats();

    return 0;
  }

  /**
   * @brief Simulate the given model, producing output at the given output times.
   *
   * Implements a simple time stepping algorithm and facilitates output.
   * 1. Calls \ref timeStep() to move the simulation forward by \ref step_size_ (starting at `params.starting_step`).
   * 2. Consults the `step_controller` to see if the step is accepted, and what the next \ref step_size_ should be.
   * 3. If accepted, advance simulation by adding \ref step_size_ to \ref current_time_ using Kahan summation and shifting
   * \ref y_new_ -> \ref y_cur_ -> \ref y_prev_.
   * 4. Any time we step over an output time, calculate dense coefficients if available.
   * 5. For each output time stepped over in the last step, interpolate the output at that time and call `out_cb`.
   *
   * Time is advanced to at least the final output time. The simulation may step over the final time, so if you
   * wish to restart simulation from the final output time, make sure to re-initialize the model's state and call
   * \ref initializeSimulation().
   *
   * @note This method can fail.
   *
   * @pre Must have called \ref initializeSimulation().
   *
   * @todo A higher-order Hermite interpolation can be used as a fallback for interpolation when dense coefficients don't exist.
   * Each step needs to sample the residual function at the beginning (and therefore end) of the step, which contains derivative
   * information. This derivative information can be re-used for Hermite interpolation.
   *
   * @todo It doesn't really make sense to have the error estimator separate from the step controller, since your choice of one will
   * affect your choice of the other. The error estimator should probably be inside the step controller, then accessed if needed.
   *
   * @todo Return an error when max steps is hit.
   *
   * @todo Check if current time is close enough to output time, and skip interpolation
   *
   * @todo Configure upper bound for skip lu step size increase
   *
   * @param out_times The times at which output is wanted. The simulation will stop once the final output time has been reached.
   * @param step_controller The step size controller to use during the simulation.
   * @param params The parameters to use during the simulation.
   * @param out_cb An optional function which, if provided, will be called once for each time in `out_times`. The simulation time of the
   * output is passed as the only argument. Before being called, the model will be updated with the output state and can be queried
   * separately by the callback if needed.
   * @param step_cb An optional function which, if provided, will be called once for each step the integrator takes. Information about the
   * step which was just taken is provided as the argument. Before being called, the model will be updated with the current value of the state
   * an can be queried separately by the callback if needed. Useful for debugging simulations.
   * @return An error code, with 0 as success.
   */
  template <class ScalarT, typename IdxT>
  int Rosenbrock<ScalarT, IdxT>::integrate(const std::vector<RealT>&                           out_times,
                                           StepController<RealT>&                              step_controller,
                                           Parameters                                          params,
                                           std::optional<std::function<void(RealT)>>           out_cb,
                                           std::optional<std::function<void(const StepInfo&)>> step_cb)
  {
    constexpr RealT ONE  = GridKit::ONE<RealT>;
    constexpr RealT ZERO = GridKit::ZERO<RealT>;

    skip_lu_ = false;
    skip_f_  = false;

    bool prev_accept = true;
    step_size_       = params.starting_step_;

    double next_step_size;

    // Kahan summation time buffer. The "leftover" time that was lost when trying to add h at some point that needs to be added
    // later
    double time_buffer = 0;

    // Generate output for each output time
    for (double out_time : out_times)
    {
      while (current_time_ < out_time && stats_.num_steps_ < params.max_steps_)
      {
        BUBBLE_FAIL(timeStep(current_time_, step_size_));

        double err = 0;

        if (step_controller.usesError())
        {
          if (err_norm_ == nullptr)
          {
            std::cerr << "The provided step controller requires the use of an error norm, but none was provided!\n";

            return -1;
          }

          State& err_vec = errorEstimate();
          err            = err_norm_->errorNorm(err_vec, *y_new_, *y_cur_, vector_handler_, memspace_);
        }

        StepControl<RealT> next_step = step_controller.nextStep(err,
                                                                StepControl<RealT>{
                                                                    .accept_    = prev_accept,
                                                                    .step_size_ = step_size_,
                                                                },
                                                                tab_.order_);
        prev_accept                  = next_step.accept_;
        next_step_size               = next_step.step_size_;

        if (prev_accept)
        {
          // Try to add the leftover time that we've stored up
          double step_size_adj = step_size_ + time_buffer;
          double next_time     = current_time_ + step_size_adj;

          // Kahan summation - keep track of how much of step_size_adj we weren't able to add to current_time
          // due to lack of precision
          time_buffer   = step_size_adj - (next_time - current_time_);
          current_time_ = next_time;

          // Since time is advancing, we need to re-evaluate the residual function and dense coefficients
          skip_f_                   = false;
          dense_coefficients_valid_ = false;

          stats_.num_steps_++;
          if (skip_lu_)
          {
            stats_.skip_lu_steps_.push_back(StepInfo{
                .sim_time_       = current_time_,
                .step_size_      = step_size_,
                .next_step_size_ = next_step_size,
                .err_est_        = err,
                .step_no_        = stats_.num_steps_,
                .skip_lu_        = skip_lu_,
                .skip_f_         = skip_f_,
                .accepted_       = prev_accept,
            });
          }
          stats_.min_step_ = std::min(stats_.min_step_, step_size_);
          stats_.max_step_ = std::max(stats_.max_step_, step_size_);

          // Shift y_new_ -> y_cur_ -> y_prev_
          // Then y_new_ is free to be replaced (contains old y_prev_)
          std::swap(y_prev_, y_cur_);
          std::swap(y_cur_, y_new_);
        }
        else
        {
          skip_f_ = true;

          stats_.rejections_.push_back(StepInfo{
              .sim_time_       = current_time_,
              .step_size_      = step_size_,
              .next_step_size_ = next_step_size,
              .err_est_        = err,
              .step_no_        = stats_.num_steps_,
              .skip_lu_        = skip_lu_,
              .skip_f_         = skip_f_,
              .accepted_       = prev_accept,
          });
        }

        // Check if we can use time delay Jacobians. If we would have increased the step size (but not too much),
        // instead keep the step size the same and use time-delay Jacobian.
        // TODO: configure upper bound here
        double step_gain = next_step_size / step_size_;
        if (params.skip_lu_ && step_gain >= ONE && step_gain <= 1.2)
        {
          skip_lu_ = true;
        }
        else
        {
          skip_lu_        = false;
          prev_step_size_ = step_size_;
          step_size_      = next_step_size;
        }

        // If there is a step_cb, update the model state and call it
        if (step_cb)
        {
          BUBBLE_FAIL(y_cur_->copyToExternal(model_->y().data(), memspace_, memspace_));
          model_->updateTime(current_time_, ZERO);

          (*step_cb)(StepInfo{
              .sim_time_       = current_time_,
              .step_size_      = step_size_,
              .next_step_size_ = next_step_size,
              .err_est_        = err,
              .step_no_        = stats_.num_steps_,
              .skip_lu_        = skip_lu_,
              .skip_f_         = skip_f_,
              .accepted_       = prev_accept,
          });
        }
      }

      // Check to make sure integration was paused because we reached an output time.
      // Other reasons, like hitting max step count, shouldn't generate output.
      if (current_time_ >= out_time)
      {
        // Theta = (t - t0) / h = (t - t1) / h + 1
        // current_time_ is t1 here
        RealT theta = (out_time - current_time_) / prev_step_size_ + ONE;

        // Generate output at the appropriate time.
        if (tab_.hasDenseOutput())
        {
          if (!dense_coefficients_valid_)
          {
            BUBBLE_FAIL(calcDenseCoeff());
            dense_coefficients_valid_ = true;
          }

          BUBBLE_FAIL(interpDense(theta));
        }
        else
        {
          // TODO: Put code for alternative interpolation (Abdou) here
          BUBBLE_FAIL(y_interp_->copyFromExternal(*y_prev_, memspace_, memspace_));
          vector_handler_.scal(1 - theta, y_interp_.get(), memspace_);
          vector_handler_.axpy(theta, y_cur_.get(), y_interp_.get(), memspace_);
        }

        // Update model with output and call out_cb if it exists
        if (out_cb)
        {
          BUBBLE_FAIL(y_interp_->copyToExternal(model_->y().data(), memspace_, memspace_));
          model_->updateTime(out_time, ZERO);

          (*out_cb)(out_time);
        }
      }
      else
      {
        BUBBLE_FAIL(y_interp_->copyFromExternal(*y_cur_, memspace_, memspace_));
        break;
      }
    }

    return 0;
  }

  /**
   * @brief Advance the simulation forward by one step, storing the new state in \ref y_new_.
   *
   * Apply the Rosenbrock scheme using the stored tableau. Each stage \f(u_i\f) is calculated as
   *
   * \f[\left(J - \frac{1}{h\gamma} M\right)u_i = -f \left(t_0 + \alpha_ih, y_0 + \sum_{j = 1}^{i - 1} a_{ij}u_j\right) - M \sum_{j=1}^{i-1} \left(\frac{c_{ij}}{h}\right)u_j,\f]
   *
   * and the next state \f(y_1\f) is calculated as
   *
   * \f[y_1 = y_0 + \sum_{j = 1}^sm_ju_j\f]
   *
   * where the coefficients \f(\gamma, \alpha_i, a_{ij}, c_{ij}, m_j\f) come from the tableau, and the mass matrix \f(M\f) and Jacobian \f(J\f) come from the model.
   *
   * This method uses some state which is maintained between calls for future calls to `timeStep()` and to communicate with \ref integrate():
   * - \ref y_cur_ is used as \f(y_0\f) and \ref y_new_ is used as \f(y_1\f)
   * - \ref asum_ is used as \f(y_0 + \sum_{j = 1}^{i - 1} a_{ij}u_j\f) for every stage except the first. Often, \f(a_{ij} = a_{i-1,j}\f), so this variable
   * can be re-used between stages to save computation. Similarly, often \f(a_{ij} = m_j\f), so this variable can be re-used for computing \ref y_new_.
   * - \ref csum_ is used as \f(M \sum_{j=1}^{i-1} \left(\frac{c_{ij}}{h}\right)u_j\f) for every stage except the first
   * - \ref RHS_ is used as \f(-f \left(t_0 + \alpha_ih, y_0 + \sum_{j = 1}^{i - 1} a_{ij}u_j\right) - M \sum_{j=1}^{i-1} \left(\frac{c_{ij}}{h}\right)u_j\f), i.e.
   * the residual evaluated at \ref asum_ plus \ref csum_, for every stage but the first. This is used as the right-hand side vector for the linear solve to solve
   * for the stage value \f(u_i\f).
   * - \ref RHS_first_stage_ is used as a special \ref RHS_ for the first stage, since it may need to be saved for the next step for the \ref skip_f_ flag.
   * Since \f(\alpha_1 = 0\f) always, this will have a value of \f(f \left(t_0, y_0\right)\f).
   * - \ref stages_ stores all of the stages \f(u_i\f). These stages are necessary to be used in future calls to \ref errorEstimate() and \ref calcDenseCoeff().
   * - \ref skip_lu_ is used as a flag set by \ref integrate() to indicate that it is appropriate to use a time-delay Jacobian by re-using the factorization
   * of the last step.
   * - \ref skip_f_ is used as a flag set by \ref integrate() to indicate that \f(t_0, y_0\f) have not changed since the last time `timeStep()` was called,
   * so \ref RHS_first_stage_ can be re-used from the previous step. This will only happen when a step was rejected, so \ref skip_lu_ should always be false,
   * and the entire first stage can't be re-used.
   * - \ref jacobian_analyzed_ keeps track of whether the Jacobian has been factored in a previous call to `timeStep()`. If so, then it can be re-factored
   * in a faster way. The first factor must be done on actual data, so it cannot be performed pre-simulation.
   *
   * @pre Must call \ref initializeSimulation() beforehand.
   *
   * @note This method can fail.
   *
   * @todo This currently does not work with non-autonomous models. Some thought needs to be put in for how we want non-autonomous models to work.
   *
   * @todo It doesn't really make sense to pass `t0` here, which may be inconsistent with \ref current_time_.
   * This function should just use \ref current_time_ in place of `t0`.
   *
   * @todo It doesn't really make sense for this method to be `public` since it relies on proper setup via \ref integrate().
   *
   * @todo May be able to move the copying of model jacobian pointers to \ref allocate().
   *
   * @todo Since \ref csum_ is multiplied by the mass matrix (which is currently diagonal with 0s for algebraic components), we can save some computation
   * by only computing the differential parts of \ref csum_ and storing it in the same variable as \ref RHS_.
   *
   * @param t0 \f(t_0\f) in the above formula
   * @param dt \f(h\f) in the above formula. The next state \f(y_1\f) will be an estimate of the state at \f(t_1 = t_0 + h\f).
   * @return An error code, with 0 as success.
   */
  template <class ScalarT, typename IdxT>
  int Rosenbrock<ScalarT, IdxT>::timeStep(RealT t0, RealT dt)
  {
    constexpr RealT ZERO      = GridKit::ZERO<RealT>;
    constexpr RealT MINUS_ONE = GridKit::MINUS_ONE<RealT>;

    // A flag to keep track of if y0 (stored in y_cur_) has been copied in to the model already, to avoid double-copying
    // for evaluating the Jacobian and residual on stage 1 (both evaluated at y0).
    bool y0_copied = false;

    // Form the left-hand side of the system. This is constant between stages.
    // Can sometimes be skipped if the method allows for time-delay Jacobians (such as w-methods).
    [[likely]]
    if (!tab_.is_w_ || !skip_lu_)
    {
      BUBBLE_FAIL(y_cur_->copyToExternal(model_->y().data(), memspace_, memspace_));
      y0_copied = true;

      // GridKit, like IDA, expects to evaluate the Jacobian J = df/dy + alpha * df/dy',
      // so we need a negative here since df/dy' = M.
      model_->updateTime(t0, MINUS_ONE / (dt * tab_.gamma_));
      BUBBLE_FAIL(model_->evaluateJacobian());
      GridKit::LinearAlgebra::CsrMatrix<RealT, IdxT>* model_jacobian = model_->getCsrJacobian();

      // TODO: This can likely be moved to allocate? These pointers should be consistent throughout the simulation
      BUBBLE_FAIL(workspace_.jacobian_->setDataPointers(
          model_jacobian->getRowData(),
          model_jacobian->getColData(),
          model_jacobian->getValues(),
          GridKit::memory::memorySpaceAsResolve(memspace_)));

      // We must factorize first (slower) and then can re-factorize (faster) on later steps
      [[likely]]
      if (workspace_.jacobian_analyzed_)
      {
        BUBBLE_FAIL(lin_solver_.refactorize());
      }
      else
      {
        BUBBLE_FAIL(lin_solver_.factorize());
        workspace_.jacobian_analyzed_ = true;
      }

      stats_.jac_evals_++;
    }

    // First stage
    [[unlikely]]
    if (skip_f_)
    {
      stats_.f_skipped_++;
    }
    else
    {
      // TODO: non-autonomous model
      if (!y0_copied)
      {
        BUBBLE_FAIL(y_cur_->copyToExternal(model_->y().data(), memspace_, memspace_));
        y0_copied = true;
      }
      model_->updateTime(t0, ZERO);
      BUBBLE_FAIL(model_->evaluateResidual());
      BUBBLE_FAIL(workspace_.RHS_first_stage_->copyFromExternal(model_->getResidual().data(), memspace_, memspace_));
      vector_handler_.scal(-1, workspace_.RHS_first_stage_.get(), memspace_);

      stats_.f_evals_++;
    }
    resolve_rhs_->setData(workspace_.RHS_first_stage_->getData(), GridKit::memory::memorySpaceAsResolve(memspace_));
    resolve_lhs_->setData(workspace_.stages_[0]->getData(), GridKit::memory::memorySpaceAsResolve(memspace_));
    BUBBLE_FAIL(lin_solver_.solve(resolve_rhs_.get(), resolve_lhs_.get()));
    stats_.decomp_solves_++;

    // Rest of stages
    for (size_t i = 1; i < tab_.num_stages_; i++)
    {
      // Calculate asum
      // We can sometimes reuse asum from the previous stage
      if (i > 1 && tab_.canReuseAsum(i))
      {
        if (tab_.A_[tab_.num_stages_ * i + i - 1] != ZERO)
          vector_handler_.axpy(tab_.A_[tab_.num_stages_ * i + i - 1], workspace_.stages_[i - 1].get(), workspace_.asum_.get(), memspace_);
      }
      else
      {
        BUBBLE_FAIL(workspace_.asum_->copyFromExternal(*y_cur_, memspace_, memspace_));

        for (size_t j = 0; j < i; j++)
        {
          if (tab_.A_[tab_.num_stages_ * i + j] != ZERO)
            vector_handler_.axpy(tab_.A_[tab_.num_stages_ * i + j], workspace_.stages_[j].get(), workspace_.asum_.get(), memspace_);
        }
      }

      // Calculate csum
      // TODO: Since csum is multiplied by the mass matrix, we can reduce calculations by just not calculating some indices
      BUBBLE_FAIL(workspace_.csum_->setToZero(memspace_));
      for (size_t j = 0; j < i; j++)
      {
        if (tab_.C_[i * tab_.num_stages_ + j] != ZERO)
        {
          vector_handler_.axpy(tab_.C_[i * tab_.num_stages_ + j] / dt, workspace_.stages_[j].get(), workspace_.csum_.get(), memspace_);
        }
      }

      // TODO: non-autonomous model
      BUBBLE_FAIL(workspace_.asum_->copyToExternal(model_->y().data(), memspace_, memspace_));
      model_->updateTime(t0 + tab_.alpha_sum_[i] * dt, ZERO);
      BUBBLE_FAIL(model_->evaluateResidual());
      workspace_.RHS_->copyFromExternal(model_->getResidual().data(), memspace_, memspace_);

      vector_handler_.scal(MINUS_ONE, workspace_.RHS_.get(), memspace_);
      vector_handler_.scal(workspace_.mass_.get(), workspace_.csum_.get(), memspace_);
      vector_handler_.axpy(MINUS_ONE, workspace_.csum_.get(), workspace_.RHS_.get(), memspace_);

      resolve_rhs_->setData(workspace_.RHS_->getData(), GridKit::memory::memorySpaceAsResolve(memspace_));
      resolve_lhs_->setData(workspace_.stages_[i]->getData(), GridKit::memory::memorySpaceAsResolve(memspace_));
      BUBBLE_FAIL(lin_solver_.solve(resolve_rhs_.get(), resolve_lhs_.get()));
      stats_.f_evals_++;
      stats_.decomp_solves_++;
    }

    // Compute the solution at time t + dt
    // It happens often where the solution is just asum of the last stage
    // plus some multiple of the last stage. In that case we can avoid a matmul
    if (tab_.canReuseAsumForOut())
    {
      std::swap(workspace_.asum_, y_new_);
      vector_handler_.axpy(tab_.m_[tab_.num_stages_ - 1], workspace_.stages_[tab_.num_stages_ - 1].get(), y_new_.get(), memspace_);
    }
    else
    {
      BUBBLE_FAIL(y_new_->copyFromExternal(*y_cur_, memspace_, memspace_));

      for (size_t j = 0; j < tab_.num_stages_; j++)
      {
        if (tab_.m_[j] != ZERO)
        {
          vector_handler_.axpy(tab_.m_[j], workspace_.stages_[j].get(), y_new_.get(), memspace_);
        }
      }
    }

    return 0;
  }

  /**
   * @brief Calculates an estimation of the error produced by the last call to \ref timeStep() which can be used as the
   * `err` argument for \ref ErrorNorm::errorNorm().
   *
   * Calculate the embedded error as
   *
   * \f[\hat{e} = \sum_{j = 1}^s e_j u_j,\f]
   *
   * where \f(e_j\f) are tableau coefficients and \f(u_j\f) are stages computed by \ref timeStep(). It happens often that
   * \f(e_j = 0\f) for all but one stage, where \f(e_j = 1\f). In that case, the stage itself is used as the error estimate,
   * and extra calculation can be avoided. For this reason, a reference to the estimate is returned to avoid an unnecessary copy.
   *
   * @pre Must call \ref timeStep(), and tableau must have coefficients for an embedded error estimator.
   * @note This method can fail.
   *
   * @todo This function is fallible, but the return type makes it difficult to return an error code. Right now it will throw an
   * error, but it should be refactored to allow returning an error code, such as by returning `std::variant`.
   *
   * @return A reference to the estimated error.
   */
  template <class ScalarT, typename IdxT>
  Rosenbrock<ScalarT, IdxT>::State& Rosenbrock<ScalarT, IdxT>::errorEstimate() const
  {
    // Test to see if the tableau allows us to use a stage as the error estimate,
    // avoiding extra computation.
    std::optional<size_t> err_stage = tab_.errorEstimatorStage();

    if (err_stage)
    {
      return *workspace_.stages_[*err_stage];
    }
    else
    {
      // TODO: could make this function return recoverable errors by using std::variant
      int err_code = workspace_.err_est_->copyFromExternal(*workspace_.stages_[0], memspace_, memspace_);

      if (err_code)
      {
        throw std::format("ReSolve::vector::Vector::copyFromExternal failed with error code {}", err_code);
      }

      vector_handler_.scal(tab_.e_[0], workspace_.err_est_.get(), memspace_);
      for (size_t j = 1; j < tab_.num_stages_; j++)
      {
        if (tab_.e_[j] != GridKit::ZERO<RealT>)
        {
          vector_handler_.axpy(tab_.e_[j], workspace_.stages_[j].get(), workspace_.err_est_.get(), memspace_);
        }
      }

      return *workspace_.err_est_;
    }
  }

  /**
   * @brief Calculate the interpolation nodes used by \ref interpDense() based on \ref stages_ computed by
   * the last call to \ref timeStep(). Nodes are stored in \ref dense_coeff_. Only needs to be called once,
   * but can be invalidated by future calls to \ref timeStep(). \ref integrate() keeps track of \ref dense_coefficients_valid_,
   * which tells you if this function is needed to be called.
   *
   * @pre Must call \ref timeStep() first.
   *
   * @note This method can fail.
   *
   * @return An error code, with 0 as success.
   */
  template <class ScalarT, typename IdxT>
  int Rosenbrock<ScalarT, IdxT>::calcDenseCoeff()
  {
    if (tab_.order_ > 2)
    {
      for (size_t j = 0; j < static_cast<size_t>(tab_.order_ - 2); j++)
      {
        BUBBLE_FAIL(dense_coeff_[j]->setToZero(memspace_));
      }

      for (size_t i = 0; i < tab_.num_stages_; i++)
      {
        for (size_t j = 0; j < static_cast<size_t>(tab_.order_ - 2); j++)
        {
          vector_handler_.axpy(tab_.H_[j * tab_.num_stages_ + i], workspace_.stages_[i].get(), dense_coeff_[j].get(), memspace_);
        }
      }
    }

    return 0;
  }

  /**
   * @brief Calculate an interpolated state at \f(\theta = \frac{t - t_0}{h}\f) between the initial state
   * \f(y_0\f) and final state \f(y_1\f) of the last step taken using dense interpolation nodes calculated
   * by \ref calcDenseCoeff().
   *
   * For a valid interpolation of appropriate order, \f(\theta\f) must be in \f([0, 1]\f), although values
   * beyond 1 can be used for extrapolation if desired.
   *
   * Uses a number of interpolation nodes equal to the order. Since \f(y_0\f) and \f(y_1\f) are interpolation
   * nodes, this method will only access \ref dense_coeff_ if the method's order is greater than 2.
   *
   * The inteporlation is calculated as
   *
   * \f[y(\theta) = (1 - \theta) y_0 + \theta \left(y_1 + (1 - \theta) \sum_{i = 1}^{p-2} \theta^{i-1} \hat{y}_i\right),\f]
   *
   * where \f(\hat{y}_i\f) are the dense interpolation nodes calculated in \ref calcDenseCoeff() and \f(p\f) is the order of the method.
   * This calculation is carried out using synthetic division.
   *
   * @pre Must call \ref calcDenseCoeff() first.
   *
   * @note This method can fail.
   *
   * @note If `theta` is very close to 0 or very close to 1, it might be better to simply use \ref y_prev_ or \ref y_cur_ instead of
   * calculating this interpolation.
   *
   * @param theta The fraction of time during the last step taken to calculate the interpolation at. \f(\theta = \frac{t - t_0}{h}\f)
   * @return An error code, with 0 as success.
   */
  template <class ScalarT, typename IdxT>
  int Rosenbrock<ScalarT, IdxT>::interpDense(RealT theta)
  {
    constexpr RealT ONE = GridKit::ONE<RealT>;

    if (tab_.order_ > 2)
    {
      BUBBLE_FAIL(y_interp_->copyFromExternal(*dense_coeff_[tab_.order_ - 3], memspace_, memspace_));

      for (size_t i = 1; i < static_cast<size_t>(tab_.order_ - 2); i++)
      {
        vector_handler_.scal(theta, y_interp_.get(), memspace_);
        vector_handler_.axpy(ONE, dense_coeff_[tab_.order_ - 3 - i].get(), y_interp_.get(), memspace_);
      }

      // TODO: This scal can be removed and absorbed into the next axpy, except that it currently isn't possible to put a scalar
      // multiple on the y term.
      vector_handler_.scal(ONE - theta, y_interp_.get(), memspace_);
    }
    else
    {
      BUBBLE_FAIL(y_interp_->setToZero(memspace_));
    }

    vector_handler_.axpy(ONE, y_cur_.get(), y_interp_.get(), memspace_);
    vector_handler_.scal(theta, y_interp_.get(), memspace_);
    vector_handler_.axpy(ONE - theta, y_prev_.get(), y_interp_.get(), memspace_);

    return 0;
  }

  /**
   * @brief Standard textbook adaptive controller. Accept if `err <= 1` and use
   *
   * \f[h_{new} = h * \min \left\{fac_{max}, \max\left\{fac_{min}, fac_{scale} \cdot e ^{-1/p}\right\}\right\}.\f]
   *
   */
  template <typename RealT>
  StepControl<RealT> AdaptiveStep<RealT>::nextStep(RealT err, StepControl<RealT> prev_step, uint8_t method_order)
  {
    StepControl<RealT> next_step = prev_step;

    double h_mult = std::min(params_.fac_max_, std::max(params_.fac_scale_ * std::pow(err, GridKit::MINUS_ONE<RealT> / method_order), params_.fac_min_));

    next_step.accept_     = err <= 1;
    next_step.step_size_ *= h_mult;

    return next_step;
  }

  /**
   * @brief Fixed step - accept every step, no matter the error, and keep the step size the same.
   *
   */
  template <typename RealT>
  StepControl<RealT> FixedStep<RealT>::nextStep([[maybe_unused]] RealT err, StepControl<RealT> prev_step, [[maybe_unused]] uint8_t method_order)
  {
    return StepControl<RealT>{
        .accept_    = true,
        .step_size_ = prev_step.step_size_,
    };
  }

  /**
   * @brief Calculate the infinity error norm as
   *
   * \f[e = \max_i\frac{|\hat{e}_i|}{Atol_i + Rtol \cdot \max\{|y_{0i}|, |y_{1i}|\}},\f]
   *
   * where \f(y_0\f) is the initial state, \f(y_1\f) is the next state, and \f(\hat{e}\f) is the estimated error made in calculating
   * the next state (typically \f(\hat{e} = y_1 - \hat{y}_1\f) for some different-order approximation \f(\hat{y}_1\f)).
   *
   * @param err \f(\hat{e}\f) in the above formula.
   * @param y \f(y_1\f) in the above formula.
   * @param yprev \f(y_0\f) in the above formula.
   * @param handler The handler to be used for performing linear algebra operations.
   * @param memspace The memory space to be used for performing linear lagebra operations.
   * @see `Rosenbrock::errorEstimate()`
   */
  template <class ScalarT, typename IdxT>
  InfNorm<ScalarT, IdxT>::RealT InfNorm<ScalarT, IdxT>::errorNorm(State& err, State& y, State& yprev, GridKit::LinearAlgebra::VectorHandler<ScalarT, IdxT>& handler, ReSolve::memory::MemorySpace memspace) const
  {
    if (int err_code = workspace_.out_->copyFromExternal(&err, memspace, memspace))
    {
      throw std::format("ReSolve::vector::Vector::copyFromExternal failed with error code {}", err_code);
    }
    if (int err_code = workspace_.scale_->copyFromExternal(&y, memspace, memspace))
    {
      throw std::format("ReSolve::vector::Vector::copyFromExternal failed with error code {}", err_code);
    }
    if (int err_code = workspace_.yprev_abs_->copyFromExternal(&yprev, memspace, memspace))
    {
      throw std::format("ReSolve::vector::Vector::copyFromExternal failed with error code {}", err_code);
    }

    handler.abs(workspace_.scale_.get(), workspace_.scale_.get(), memspace);
    handler.abs(workspace_.yprev_abs_.get(), workspace_.scale_.get(), memspace);
    handler.max(workspace_.yprev_abs_.get(), workspace_.scale_.get(), workspace_.scale_.get(), memspace);

    // TODO: This scal shouldn't be necessary, but axpy doesn't support scaling the y parameter. In the future,
    // the scaling should be able to be put on the next axpy.
    handler.scal(params_.rel_tol_, workspace_.scale_.get(), memspace);
    handler.axpy(GridKit::ONE<RealT>, params_.abs_tol_.get(), workspace_.scale_.get(), memspace);
    handler.diagSolve(workspace_.scale_.get(), workspace_.out_.get(), memspace);

    return handler.amax(workspace_.out_.get(), memspace);
  }

  template class Rosenbrock<double, int>;

  template class FixedStep<double>;

  template class AdaptiveStep<double>;
} // namespace Integrator
