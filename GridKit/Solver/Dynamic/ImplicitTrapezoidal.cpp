#include "ImplicitTrapezoidal.hpp"

#include <algorithm>
#include <cmath>

#include <GridKit/Constants.hpp>

#define BUBBLE_FAIL(arg) \
  do                     \
  {                      \
    if (int err = (arg)) \
    {                    \
      return err;        \
    }                    \
  } while (false)

#define BUBBLE_STEP_FAIL(arg)        \
  do                                 \
  {                                  \
    if (int err = (arg))             \
    {                                \
      jacobian_valid_ = false;       \
      (void) restoreAcceptedState(); \
      return err;                    \
    }                                \
  } while (false)

namespace AnalysisManager
{
  namespace NativeDynamicSolver
  {
    template <class ScalarT, typename IdxT>
    ImplicitTrapezoidal<ScalarT, IdxT>::ImplicitTrapezoidal(
        GridKit::Model::Evaluator<ScalarT, IdxT>*             model,
        GridKit::LinearAlgebra::LinearSolver<ScalarT, IdxT>&  linear_solver,
        GridKit::LinearAlgebra::VectorHandler<ScalarT, IdxT>& vector_handler,
        const ErrorNorm<ScalarT, IdxT>*                       error_norm,
        GridKit::memory::MemorySpace                          memspace)
      : model_(model), linear_solver_(linear_solver), vector_handler_(vector_handler), error_norm_(error_norm), memspace_(memspace)
    {
    }

    template <class ScalarT, typename IdxT>
    int ImplicitTrapezoidal<ScalarT, IdxT>::allocate()
    {
      const size_t size = static_cast<size_t>(model_->size());

      current_state_       = std::make_unique<State>(size);
      current_derivative_  = std::make_unique<State>(size);
      previous_state_      = std::make_unique<State>(size);
      previous_derivative_ = std::make_unique<State>(size);
      iterate_             = std::make_unique<State>(size);
      residual_            = std::make_unique<State>(size);
      correction_          = std::make_unique<State>(size);

      BUBBLE_FAIL(current_state_->allocate(memspace_));
      BUBBLE_FAIL(current_derivative_->allocate(memspace_));
      BUBBLE_FAIL(previous_state_->allocate(memspace_));
      BUBBLE_FAIL(previous_derivative_->allocate(memspace_));
      BUBBLE_FAIL(iterate_->allocate(memspace_));
      BUBBLE_FAIL(residual_->allocate(memspace_));
      BUBBLE_FAIL(correction_->allocate(memspace_));
      BUBBLE_FAIL(correction_->setToZero(memspace_));

      return 0;
    }

    template <class ScalarT, typename IdxT>
    int ImplicitTrapezoidal<ScalarT, IdxT>::initializeSimulation(RealT t0)
    {
      if (error_norm_ == nullptr || !model_->hasJacobian() || model_->getCsrJacobian() == nullptr)
      {
        return -1;
      }

      BUBBLE_FAIL(current_state_->copyFromExternal(model_->y(), memspace_, memspace_));
      BUBBLE_FAIL(current_derivative_->copyFromExternal(model_->yp(), memspace_, memspace_));
      BUBBLE_FAIL(linear_solver_.configureSolver(*model_->getCsrJacobian()));

      current_time_        = t0;
      previous_time_       = t0;
      jacobian_factorized_ = false;
      jacobian_valid_      = false;
      jacobian_step_size_  = GridKit::ZERO<RealT>;
      stats_               = Stats{};
      initialized_         = true;
      model_->updateTime(t0, GridKit::ZERO<RealT>);
      return 0;
    }

    template <class ScalarT, typename IdxT>
    bool ImplicitTrapezoidal<ScalarT, IdxT>::validParameters(const Parameters& p) const
    {
      return p.step_size_ > 0
             && p.max_steps_ > 0
             && p.max_newton_iterations_ > 0
             && p.armijo_constant_ > 0
             && p.armijo_constant_ < 1
             && p.backtrack_factor_ > 0
             && p.backtrack_factor_ < 1
             && p.minimum_step_length_ > 0
             && p.minimum_step_length_ <= 1
             && p.stagnation_ratio_ > 0
             && p.stagnation_ratio_ < 1
             && p.max_stagnant_iterations_ > 0
             && p.max_jacobian_refreshes_ > 0;
    }

    template <class ScalarT, typename IdxT>
    int ImplicitTrapezoidal<ScalarT, IdxT>::restoreAcceptedState()
    {
      BUBBLE_FAIL(model_->y().copyFromExternal(*current_state_, memspace_, memspace_));
      BUBBLE_FAIL(model_->yp().copyFromExternal(*current_derivative_, memspace_, memspace_));
      model_->updateTime(current_time_, GridKit::ZERO<RealT>);
      return 0;
    }

    template <class ScalarT, typename IdxT>
    int ImplicitTrapezoidal<ScalarT, IdxT>::evaluateStepResidual(RealT t1, RealT dt, State& residual)
    {
      BUBBLE_FAIL(model_->yp().copyFromExternal(*current_derivative_, memspace_, memspace_));
      vector_handler_.scal(-1, &model_->yp(), memspace_);
      vector_handler_.axpy(-2 / dt, current_state_.get(), &model_->yp(), memspace_);
      vector_handler_.axpy(2 / dt, &model_->y(), &model_->yp(), memspace_);

      model_->updateTime(t1, GridKit::ZERO<RealT>);
      BUBBLE_FAIL(model_->evaluateResidual());
      stats_.num_residual_evaluations_++;
      return residual.copyFromExternal(model_->getResidual(), memspace_, memspace_);
    }

    template <class ScalarT, typename IdxT>
    int ImplicitTrapezoidal<ScalarT, IdxT>::evaluateStepJacobian(RealT t1, RealT dt)
    {
      model_->updateTime(t1, RealT(2) / dt);
      BUBBLE_FAIL(model_->evaluateJacobian());
      stats_.num_jacobian_evaluations_++;
      BUBBLE_FAIL(linear_solver_.setupSolver(jacobian_factorized_));
      jacobian_factorized_ = true;
      jacobian_valid_      = true;
      jacobian_step_size_  = dt;
      return 0;
    }

    template <class ScalarT, typename IdxT>
    typename ImplicitTrapezoidal<ScalarT, IdxT>::RealT
    ImplicitTrapezoidal<ScalarT, IdxT>::residualMerit(State& residual) const
    {
      return RealT(0.5) * vector_handler_.dot(&residual, &residual, memspace_);
    }

    template <class ScalarT, typename IdxT>
    int ImplicitTrapezoidal<ScalarT, IdxT>::timeStep(RealT t0, RealT dt, const Parameters& parameters)
    {
      if (!initialized_ || !validParameters(parameters) || dt <= 0 || t0 != current_time_)
      {
        return -1;
      }

      const RealT t1 = t0 + dt;
      BUBBLE_STEP_FAIL(iterate_->copyFromExternal(*current_state_, memspace_, memspace_));
      vector_handler_.axpy(ScalarT(dt), current_derivative_.get(), iterate_.get(), memspace_);
      BUBBLE_STEP_FAIL(model_->y().copyFromExternal(*iterate_, memspace_, memspace_));
      BUBBLE_STEP_FAIL(evaluateStepResidual(t1, dt, *residual_));

      RealT  merit               = residualMerit(*residual_);
      bool   converged           = merit == 0;
      size_t step_iterations     = 0;
      size_t step_backtracks     = 0;
      size_t step_refreshes      = 0;
      size_t stagnant_iterations = 0;

      auto refreshJacobian = [&]() -> int
      {
        BUBBLE_FAIL(evaluateStepJacobian(t1, dt));
        step_refreshes++;
        stats_.num_jacobian_refreshes_++;
        stagnant_iterations = 0;
        return 0;
      };

      if (!converged)
      {
        if (jacobian_valid_ && jacobian_step_size_ == dt)
        {
          stats_.num_jacobian_reuses_++;
        }
        else
        {
          BUBBLE_STEP_FAIL(evaluateStepJacobian(t1, dt));
        }
      }

      while (!converged && step_iterations < parameters.max_newton_iterations_)
      {
        step_iterations++;
        const RealT old_merit = merit;
        vector_handler_.scal(ScalarT(-1), residual_.get(), memspace_);
        BUBBLE_STEP_FAIL(linear_solver_.solve(*residual_, *correction_));
        stats_.num_linear_solves_++;
        stats_.num_newton_iterations_++;

        if (error_norm_->errorNorm(*correction_, model_->y(), *current_state_, vector_handler_, memspace_) <= 1)
        {
          converged = true;
          break;
        }

        RealT length   = 1;
        bool  accepted = false;
        while (length >= parameters.minimum_step_length_)
        {
          BUBBLE_STEP_FAIL(model_->y().copyFromExternal(*iterate_, memspace_, memspace_));
          vector_handler_.axpy(length, correction_.get(), &model_->y(), memspace_);
          BUBBLE_STEP_FAIL(evaluateStepResidual(t1, dt, *residual_));
          merit = residualMerit(*residual_);
          if (merit <= (RealT(1) - parameters.armijo_constant_ * length) * old_merit)
          {
            accepted = true;
            break;
          }
          length *= parameters.backtrack_factor_;
          step_backtracks++;
          stats_.num_backtracks_++;
        }

        if (!accepted)
        {
          if (step_refreshes >= parameters.max_jacobian_refreshes_)
          {
            break;
          }

          BUBBLE_STEP_FAIL(model_->y().copyFromExternal(*iterate_, memspace_, memspace_));
          BUBBLE_STEP_FAIL(evaluateStepResidual(t1, dt, *residual_));
          merit = old_merit;
          BUBBLE_STEP_FAIL(refreshJacobian());
          continue;
        }

        vector_handler_.scal(length, correction_.get(), memspace_);
        converged = error_norm_->errorNorm(*correction_, model_->y(), *current_state_, vector_handler_, memspace_) <= 1
                    || merit == 0;
        BUBBLE_STEP_FAIL(iterate_->copyFromExternal(model_->y(), memspace_, memspace_));

        if (!converged)
        {
          stagnant_iterations = merit > parameters.stagnation_ratio_ * old_merit ? stagnant_iterations + 1 : 0;
          if (stagnant_iterations >= parameters.max_stagnant_iterations_)
          {
            if (step_refreshes >= parameters.max_jacobian_refreshes_)
            {
              break;
            }

            BUBBLE_STEP_FAIL(refreshJacobian());
          }
        }
      }

      if (!converged)
      {
        stats_.num_convergence_failures_++;
        jacobian_valid_ = false;
        BUBBLE_FAIL(restoreAcceptedState());
        return -1;
      }

      std::swap(previous_state_, current_state_);
      std::swap(previous_derivative_, current_derivative_);
      previous_time_ = current_time_;
      BUBBLE_FAIL(current_state_->copyFromExternal(*iterate_, memspace_, memspace_));
      BUBBLE_FAIL(current_derivative_->copyFromExternal(model_->yp(), memspace_, memspace_));

      current_time_ = t1;
      stats_.num_steps_++;
      model_->updateTime(current_time_, GridKit::ZERO<RealT>);
      last_step_info_ = StepInfo{current_time_, dt, stats_.num_steps_, step_iterations, step_backtracks};
      return 0;
    }

    template <class ScalarT, typename IdxT>
    int ImplicitTrapezoidal<ScalarT, IdxT>::integrate(
        const std::vector<RealT>&                           output_times,
        Parameters                                          parameters,
        std::optional<std::function<void(RealT)>>           output_callback,
        std::optional<std::function<void(const StepInfo&)>> step_callback)
    {
      if (!initialized_
          || !validParameters(parameters)
          || !std::is_sorted(output_times.begin(), output_times.end())
          || (!output_times.empty() && output_times.front() < current_time_))
      {
        return -1;
      }

      for (RealT output_time : output_times)
      {
        while (current_time_ < output_time)
        {
          if (stats_.num_steps_ >= parameters.max_steps_)
          {
            BUBBLE_FAIL(restoreAcceptedState());
            return -1;
          }
          BUBBLE_FAIL(timeStep(current_time_, parameters.step_size_, parameters));
          if (step_callback)
          {
            (*step_callback)(last_step_info_);
          }
        }

        if (!output_callback)
        {
          continue;
        }

        if (output_time == current_time_)
        {
          BUBBLE_FAIL(restoreAcceptedState());
        }
        else
        {
          const RealT theta = (output_time - previous_time_) / (current_time_ - previous_time_);
          BUBBLE_FAIL(model_->y().copyFromExternal(*previous_state_, memspace_, memspace_));
          vector_handler_.scal(1 - theta, &model_->y(), memspace_);
          vector_handler_.axpy(theta, current_state_.get(), &model_->y(), memspace_);
          BUBBLE_FAIL(model_->yp().copyFromExternal(*previous_derivative_, memspace_, memspace_));
          vector_handler_.scal(1 - theta, &model_->yp(), memspace_);
          vector_handler_.axpy(theta, current_derivative_.get(), &model_->yp(), memspace_);
          model_->updateTime(output_time, GridKit::ZERO<RealT>);
        }
        (*output_callback)(output_time);
      }
      return 0;
    }

    template class ImplicitTrapezoidal<double, int>;
  } // namespace NativeDynamicSolver
} // namespace AnalysisManager

#undef BUBBLE_FAIL
#undef BUBBLE_STEP_FAIL
