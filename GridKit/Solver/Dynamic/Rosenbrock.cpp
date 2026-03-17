#include "Rosenbrock.hpp"

#include <format>
#include <iomanip>
#include <sstream>

#include <sundials/sundials_types.h>

#define BUBBLE_FAIL(arg) \
  do                     \
  {                      \
    if (int err = (arg)) \
      return err;        \
  } while (false)

namespace Integrator
{
  template <class ScalarT, typename IdxT>
  std::string Rosenbrock<ScalarT, IdxT>::StepInfo::csv_report() const
  {
    std::stringstream out;
    out << std::scientific << std::setprecision(20)
        << sim_time << ','
        << step_size << ','
        << next_step_size << ','
        << err_est << ','
        << step_no << ','
        << skip_lu << ','
        << skip_f << ','
        << accepted;

    return out.str();
  }

  template <class ScalarT, typename IdxT>
  std::string Rosenbrock<ScalarT, IdxT>::StepInfo::report() const
  {
    std::stringstream out;
    out << std::scientific << std::setprecision(20)
        << "Simulation Time: " << sim_time << '\n'
        << "Step Size:       " << step_size << '\n'
        << "Next Step Size:  " << next_step_size << '\n'
        << "Error Estimate:  " << err_est << '\n'
        << "Step Number:     " << step_no << '\n'
        << "Skip LU:         " << skip_lu << '\n'
        << "Skip F:          " << skip_f << '\n'
        << "Accepted:        " << accepted;

    return out.str();
  }

  template <class ScalarT, typename IdxT>
  std::string Rosenbrock<ScalarT, IdxT>::Stats::report() const
  {
    std::stringstream out;
    out << "Rejections: " << rejections.size()
        << "\nSteps: " << num_steps
        << "\nSkip LU Steps: " << skip_lu_steps.size()
        << "\nMin Step: " << min_step
        << "\nMax Step: " << max_step
        << "\nRHS Function Evals: " << f_evals
        << "\nRHS Function Skipped: " << f_skipped
        << "\nJacobian Evals: " << jac_evals
        << "\nLinear Solves: " << decomp_solves;

    return out.str();
  }

  template <class ScalarT, typename IdxT>
  Rosenbrock<ScalarT, IdxT>::Stats::Stats& Rosenbrock<ScalarT, IdxT>::Stats::operator+=(const Stats& other)
  {
    rejections.insert(rejections.end(), other.rejections.begin(), other.rejections.end());
    skip_lu_steps.insert(skip_lu_steps.end(), other.skip_lu_steps.begin(), other.skip_lu_steps.end());

    num_steps     += other.num_steps;
    f_evals       += other.f_evals;
    f_skipped     += other.f_skipped;
    jac_evals     += other.jac_evals;
    decomp_solves += other.decomp_solves;

    min_step = std::min(min_step, other.min_step);
    max_step = std::max(max_step, other.max_step);

    return *this;
  }

  template <class ScalarT, typename IdxT>
  constexpr bool Rosenbrock<ScalarT, IdxT>::Tableau::can_reuse_asum(size_t stage) const
  {
    assert(stage < num_stages);

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

  template <class ScalarT, typename IdxT>
  constexpr bool Rosenbrock<ScalarT, IdxT>::Tableau::can_reuse_asum_for_out() const
  {
    if (num_stages == 1)
      return false;

    for (size_t j = 0; j < num_stages - 1; j++)
    {
      if (getA(num_stages - 1, j) != m[j])
      {
        return false;
      }
    }

    return true;
  }

  template <class ScalarT, typename IdxT>
  constexpr std::tuple<bool, size_t> Rosenbrock<ScalarT, IdxT>::Tableau::error_estimator_stage() const
  {
    std::tuple<bool, size_t> re = {false, 0};
    for (size_t j = 0; j < num_stages; j++)
    {
      if (e[j] == 1.0 && !std::get<0>(re))
      {
        re = {true, j};
      }
      else if (e[j] != 0.0)
      {
        return {false, 0};
      }
    }
    return re;
  }

  template <class ScalarT, typename IdxT>
  Rosenbrock<ScalarT, IdxT>::Rosenbrock(Tableau&&                                 tab,
                                        GridKit::Model::Evaluator<ScalarT, IdxT>* model,
                                        ReSolve::SystemSolver&                    lin_solver,
                                        ReSolve::VectorHandler&                   vector_handler,
                                        const ErrorNorm*                          err_norm,
                                        ReSolve::memory::MemorySpace              memspace)
    : tab_(std::move(tab)),
      model_(model),
      lin_solver_(lin_solver),
      vector_handler_(vector_handler),
      err_norm_(err_norm),
      memspace_(memspace)
  {
  }

  template <class ScalarT, typename IdxT>
  int Rosenbrock<ScalarT, IdxT>::allocate()
  {
    size_t size = static_cast<size_t>(model_->size());

    y_prev_          = std::make_unique<State>(size);
    y_cur_           = std::make_unique<State>(size);
    y_interp_        = std::make_unique<State>(size);
    asum_            = std::make_unique<State>(size);
    csum_            = std::make_unique<State>(size);
    RHS_             = std::make_unique<State>(size);
    RHS_first_stage_ = std::make_unique<State>(size);
    dFdt_            = std::make_unique<State>(size);
    mass_            = std::make_unique<State>(size);
    y_new_           = std::make_unique<State>(size);

    BUBBLE_FAIL(y_prev_->allocate(memspace_));
    BUBBLE_FAIL(y_cur_->allocate(memspace_));
    BUBBLE_FAIL(y_interp_->allocate(memspace_));
    BUBBLE_FAIL(asum_->allocate(memspace_));
    BUBBLE_FAIL(csum_->allocate(memspace_));
    BUBBLE_FAIL(RHS_->allocate(memspace_));
    BUBBLE_FAIL(RHS_first_stage_->allocate(memspace_));
    BUBBLE_FAIL(dFdt_->allocate(memspace_));
    BUBBLE_FAIL(mass_->allocate(memspace_));
    BUBBLE_FAIL(y_new_->allocate(memspace_));

    if (tab_.e)
    {
      auto [can_use_stage, err_stage] = tab_.error_estimator_stage();
      if (!can_use_stage)
      {
        err_est_ = std::make_unique<State>(size);
        BUBBLE_FAIL(err_est_->allocate(memspace_));
      }
    }

    stages_ = std::make_unique<std::unique_ptr<State>[]>(tab_.num_stages);
    for (size_t i = 0; i < tab_.num_stages; i++)
    {
      stages_[i] = std::make_unique<State>(size);
      BUBBLE_FAIL(stages_[i]->allocate(memspace_));
    }

    if (tab_.order > 2)
    {
      dense_coeff_ = std::make_unique<std::unique_ptr<State>[]>(tab_.order - 2);
      for (size_t i = 0; i < tab_.order - 2; i++)
      {
        dense_coeff_[i] = std::make_unique<State>(size);
        BUBBLE_FAIL(dense_coeff_[i]->allocate(memspace_));
      }
    }

    jacobian_ = std::make_unique<ReSolve::matrix::Csr>(size, size, model_->getCsrJacobian()->getNnz());

    return 0;
  }

  template <class ScalarT, typename IdxT>
  int Rosenbrock<ScalarT, IdxT>::initializeSimulation(RealT t0)
  {
    current_time_ = t0;
    BUBBLE_FAIL(y_cur_->copyFromExternal(model_->y().data(), memspace_, memspace_));
    jacobian_analyzed_ = false;

    GridKit::LinearAlgebra::CsrMatrix<RealT, IdxT>* model_jacobian = model_->getCsrJacobian();
    BUBBLE_FAIL(jacobian_->setDataPointers(
        model_jacobian->getRowData(),
        model_jacobian->getColData(),
        model_jacobian->getValues(),
        memspace_));
    BUBBLE_FAIL(lin_solver_.setMatrix(jacobian_.get()));
    BUBBLE_FAIL(lin_solver_.analyze());
    BUBBLE_FAIL(lin_solver_.preconditionerSetup());

    if (model_->tag().size() != static_cast<size_t>(model_->size()))
    {
      std::cerr << "Model tag is either unset or does not match the size of the model\n";
      return 1;
    }

    std::unique_ptr<RealT[]> mass = std::make_unique<RealT[]>(model_->tag().size());
    for (size_t i = 0; i < static_cast<size_t>(model_->size()); i++)
    {
      mass[i] = model_->tag()[i] ? 1.0 : 0.0;
    }
    BUBBLE_FAIL(mass_->copyFromExternal(mass.get(), memspace_, memspace_));

    stats_ = Stats();

    return 0;
  }

  template <class ScalarT, typename IdxT>
  int Rosenbrock<ScalarT, IdxT>::integrate(const std::vector<double>&                          out_times,
                                           StepController&                                     step_controller,
                                           Parameters                                          params,
                                           std::optional<std::function<void(double)>>          out_cb,
                                           std::optional<std::function<void(const StepInfo&)>> step_cb)
  {
    skip_lu_ = false;
    skip_f_  = false;

    bool prev_accept = true;
    step_size_       = params.starting_step;

    double next_step_size;

    // Kahan summation time buffer. The "leftover" time that was lost when trying to add h at some point that needs to be added
    // later
    double time_buffer = 0;

    for (double out_time : out_times)
    {
      while (current_time_ < out_time && stats_.num_steps < params.max_steps)
      {
        BUBBLE_FAIL(time_step(current_time_, step_size_));

        double err = 0;

        if (step_controller.usesError())
        {
          if (err_norm_ == nullptr)
          {
            std::cerr << "The provided step controller requires the use of an error norm, but none was provided!\n";

            return -1;
          }

          State& err_vec = error_estimate();
          err            = err_norm_->errorNorm(err_vec, *y_new_, *y_cur_, vector_handler_, memspace_);
        }

        StepControl next_step = step_controller.nextStep(err,
                                                         StepControl{
                                                             .accept    = prev_accept,
                                                             .step_size = step_size_,
                                                         },
                                                         tab_.order);
        prev_accept           = next_step.accept;
        next_step_size        = next_step.step_size;

        if (prev_accept)
        {
          // Try to add the leftover time that we've stored up
          double step_size_adj = step_size_ + time_buffer;
          double next_time     = current_time_ + step_size_adj;

          // Kahan summation - keep track of how much of step_size_adj we weren't able to add to current_time
          // due to lack of precision
          time_buffer   = step_size_adj - (next_time - current_time_);
          current_time_ = next_time;

          // step_cb(current_time, yout);
          skip_f_ = false;

          stats_.num_steps++;
          if (skip_lu_)
          {
            stats_.skip_lu_steps.push_back(StepInfo{
                .sim_time       = current_time_,
                .step_size      = step_size_,
                .next_step_size = next_step_size,
                .err_est        = err,
                .step_no        = stats_.num_steps,
                .skip_lu        = skip_lu_,
                .skip_f         = skip_f_,
                .accepted       = prev_accept,
            });
          }
          stats_.min_step = std::min(stats_.min_step, step_size_);
          stats_.max_step = std::max(stats_.max_step, step_size_);

          std::swap(y_prev_, y_cur_);
          std::swap(y_cur_, y_new_);
          dense_coefficients_valid = false;
        }
        else
        {
          skip_f_ = true;

          stats_.rejections.push_back(StepInfo{
              .sim_time       = current_time_,
              .step_size      = step_size_,
              .next_step_size = next_step_size,
              .err_est        = err,
              .step_no        = stats_.num_steps,
              .skip_lu        = skip_lu_,
              .skip_f         = skip_f_,
              .accepted       = prev_accept,
          });
        }

        double step_gain = next_step_size / step_size_;
        if (params.skip_lu && step_gain >= 1 && step_gain <= 1.2)
        {
          skip_lu_ = true;
        }
        else
        {
          skip_lu_        = false;
          prev_step_size_ = step_size_;
          step_size_      = next_step_size;
        }

        if (step_cb)
        {
          BUBBLE_FAIL(y_cur_->copyToExternal(model_->y().data(), memspace_, memspace_));
          model_->updateTime(current_time_, 0.0);

          (*step_cb)(StepInfo{
              .sim_time       = current_time_,
              .step_size      = step_size_,
              .next_step_size = next_step_size,
              .err_est        = err,
              .step_no        = stats_.num_steps,
              .skip_lu        = skip_lu_,
              .skip_f         = skip_f_,
              .accepted       = prev_accept,
          });
        }
      }

      if (current_time_ >= out_time)
      {
        if (tab_.hasDenseOutput())
        {
          if (!dense_coefficients_valid)
          {
            BUBBLE_FAIL(calc_dense_coeff());
            dense_coefficients_valid = true;
          }

          double theta = (out_time - current_time_) / prev_step_size_ + 1;
          BUBBLE_FAIL(interp_dense(theta));
        }
        else
        {
          // TODO: Put code for alternative interpolation (Abdou) here
          double theta = (out_time - current_time_) / prev_step_size_ + 1;
          BUBBLE_FAIL(y_interp_->copyFromExternal(y_prev_.get(), memspace_, memspace_));
          vector_handler_.scal(1 - theta, y_interp_.get(), memspace_);
          vector_handler_.axpy(theta, y_cur_.get(), y_interp_.get(), memspace_);
        }

        if (out_cb)
        {
          BUBBLE_FAIL(y_interp_->copyToExternal(model_->y().data(), memspace_, memspace_));
          model_->updateTime(out_time, 0.0);

          (*out_cb)(out_time);
        }
      }
      else
      {
        BUBBLE_FAIL(y_interp_->copyFromExternal(y_cur_.get(), memspace_, memspace_));
        break;
      }
    }

    return 0;
  }

  template <class ScalarT, typename IdxT>
  int Rosenbrock<ScalarT, IdxT>::time_step(double t0, double dt)
  {
    bool y0_copied = false;

    // Form the left-hand side of the system
    // Can sometimes be skipped if the method is a w-method
    [[likely]]
    if (!tab_.is_w || !skip_lu_)
    {
      BUBBLE_FAIL(y_cur_->copyToExternal(model_->y().data(), memspace_, memspace_));
      y0_copied = true;
      model_->updateTime(t0, -1.0 / (dt * tab_.gamma));
      BUBBLE_FAIL(model_->evaluateJacobian());
      GridKit::LinearAlgebra::CsrMatrix<RealT, IdxT>* model_jacobian = model_->getCsrJacobian();
      BUBBLE_FAIL(jacobian_->setDataPointers(
          model_jacobian->getRowData(),
          model_jacobian->getColData(),
          model_jacobian->getValues(),
          memspace_));

      [[likely]]
      if (jacobian_analyzed_)
      {
        BUBBLE_FAIL(lin_solver_.refactorize());
      }
      else
      {
        BUBBLE_FAIL(lin_solver_.factorize());
      }

      stats_.jac_evals++;
    }

    // First stage
    [[unlikely]]
    if (skip_f_)
    {
      stats_.f_skipped++;
    }
    else
    {
      // TODO: non-autonomous model
      if (!y0_copied)
      {
        BUBBLE_FAIL(y_cur_->copyToExternal(model_->y().data(), memspace_, memspace_));
        y0_copied = true;
      }
      model_->updateTime(t0, 0.0);
      BUBBLE_FAIL(model_->evaluateResidual());
      BUBBLE_FAIL(RHS_first_stage_->copyFromExternal(model_->getResidual().data(), memspace_, memspace_));
      vector_handler_.scal(-1, RHS_first_stage_.get(), memspace_);

      stats_.f_evals++;
    }
    BUBBLE_FAIL(lin_solver_.solve(RHS_first_stage_.get(), stages_[0].get()));
    stats_.decomp_solves++;

    // Rest of stages
    for (size_t i = 1; i < tab_.num_stages; i++)
    {
      // Calculate asum
      // We can sometimes reuse asum from the previous stage
      if (i > 1 && tab_.can_reuse_asum(i))
      {
        if (tab_.A[tab_.num_stages * i + i - 1] != 0.0)
          vector_handler_.axpy(tab_.A[tab_.num_stages * i + i - 1], stages_[i - 1].get(), asum_.get(), memspace_);
      }
      else
      {
        BUBBLE_FAIL(asum_->copyFromExternal(y_cur_.get(), memspace_, memspace_));

        for (size_t j = 0; j < i; j++)
        {
          if (tab_.A[tab_.num_stages * i + j] != 0.0)
            vector_handler_.axpy(tab_.A[tab_.num_stages * i + j], stages_[j].get(), asum_.get(), memspace_);
        }
      }

      // Calculate csum
      // TODO: Since csum is multiplied by the mass matrix, we can reduce calculations by just not calculating some indices
      BUBBLE_FAIL(csum_->setToZero(memspace_));
      for (size_t j = 0; j < i; j++)
      {
        if (tab_.C[i * tab_.num_stages + j] != 0.0)
        {
          vector_handler_.axpy(tab_.C[i * tab_.num_stages + j] / dt, stages_[j].get(), csum_.get(), memspace_);
        }
      }

      // TODO: non-autonomous model
      BUBBLE_FAIL(asum_->copyToExternal(model_->y().data(), memspace_, memspace_));
      model_->updateTime(t0 + tab_.alpha_sum[i] * dt, 0.0);
      BUBBLE_FAIL(model_->evaluateResidual());
      RHS_->copyFromExternal(model_->getResidual().data(), memspace_, memspace_);

      vector_handler_.scal(-1, RHS_.get(), memspace_);
      vector_handler_.scal(mass_.get(), csum_.get(), memspace_);
      vector_handler_.axpy(-1, csum_.get(), RHS_.get(), memspace_);

      BUBBLE_FAIL(lin_solver_.solve(RHS_.get(), stages_[i].get()));
      stats_.f_evals++;
      stats_.decomp_solves++;
    }

    // Compute the solution at time t + dt
    // It happens often where the solution is just asum of the last stage
    // plus some multiple of the last stage. In that case we can avoid a matmul
    if (tab_.can_reuse_asum_for_out())
    {
      std::swap(asum_, y_new_);
      vector_handler_.axpy(tab_.m[tab_.num_stages - 1], stages_[tab_.num_stages - 1].get(), y_new_.get(), memspace_);
    }
    else
    {
      BUBBLE_FAIL(y_new_->copyFromExternal(y_cur_.get(), memspace_, memspace_));

      for (size_t j = 0; j < tab_.num_stages; j++)
      {
        if (tab_.m[j] != 0.0)
        {
          vector_handler_.axpy(tab_.m[j], stages_[j].get(), y_new_.get(), memspace_);
        }
      }
    }

    return 0;
  }

  template <class ScalarT, typename IdxT>
  State& Rosenbrock<ScalarT, IdxT>::error_estimate() const
  {
    auto [can_use_stage, err_stage] = tab_.error_estimator_stage();

    if (can_use_stage)
    {
      return *stages_[err_stage];
    }
    else
    {
      // TODO: could make this function return recoverable errors by using std::variant
      int err_code = err_est_->copyFromExternal(stages_[0].get(), memspace_, memspace_);

      if (err_code)
      {
        throw std::format("ReSolve::vector::Vector::copyFromExternal failed with error code {}", err_code);
      }

      vector_handler_.scal(tab_.e[0], stages_[0].get(), memspace_);
      for (size_t j = 1; j < tab_.num_stages; j++)
      {
        if (tab_.e[j] != 0.0)
        {
          vector_handler_.axpy(tab_.e[j], stages_[j].get(), err_est_.get(), memspace_);
        }
      }

      return *err_est_;
    }
  }

  template <class ScalarT, typename IdxT>
  int Rosenbrock<ScalarT, IdxT>::calc_dense_coeff()
  {
    if (tab_.order > 2)
    {
      for (size_t j = 0; j < tab_.order - 2; j++)
      {
        BUBBLE_FAIL(dense_coeff_[j]->setToZero(memspace_));
      }

      for (size_t i = 0; i < tab_.num_stages; i++)
      {
        for (size_t j = 0; j < tab_.order - 2; j++)
        {
          vector_handler_.axpy(tab_.H[j * tab_.num_stages + i], stages_[i].get(), dense_coeff_[j].get(), memspace_);
        }
      }
    }

    return 0;
  }

  template <class ScalarT, typename IdxT>
  int Rosenbrock<ScalarT, IdxT>::interp_dense(double theta)
  {
    double one = 1.0;
    if (tab_.order > 2)
    {
      BUBBLE_FAIL(y_interp_->copyFromExternal(dense_coeff_[tab_.order - 3].get(), memspace_, memspace_));

      for (size_t i = 1; i < tab_.order - 2; i++)
      {
        vector_handler_.scal(theta, y_interp_.get(), memspace_);
        vector_handler_.axpy(one, dense_coeff_[tab_.order - 3 - i].get(), y_interp_.get(), memspace_);
      }
    }
    else
    {
      BUBBLE_FAIL(y_interp_->setToZero(memspace_));
    }

    double omt = 1 - theta;

    vector_handler_.scal(omt, y_interp_.get(), memspace_);
    vector_handler_.axpy(one, y_cur_.get(), y_interp_.get(), memspace_);
    vector_handler_.scal(theta, y_interp_.get(), memspace_);
    vector_handler_.axpy(omt, y_prev_.get(), y_interp_.get(), memspace_);

    return 0;
  }

  StepControl AdaptiveStep::nextStep(double err, StepControl prev_step, uint8_t method_order)
  {
    StepControl next_step = prev_step;

    double h_mult = std::min(params_.facmax, std::max(params_.facscale * std::pow(err, -1.0 / method_order), params_.facmin));

    next_step.accept     = err <= 1;
    next_step.step_size *= h_mult;

    return next_step;
  }

  StepControl FixedStep::nextStep([[maybe_unused]] double err, [[maybe_unused]] StepControl prev_step, [[maybe_unused]] uint8_t method_order)
  {
    return StepControl{
        .accept    = true,
        .step_size = prev_step.step_size,
    };
  }

  double InfNorm::errorNorm(State& err, State& y, State& yprev, ReSolve::VectorHandler& handler, ReSolve::memory::MemorySpace memspace) const
  {
    workspace_.out_->copyFromExternal(&err, memspace, memspace);
    workspace_.scale_->copyFromExternal(&y, memspace, memspace);
    workspace_.yprev_abs_->copyFromExternal(&yprev, memspace, memspace);

    handler.abs(workspace_.scale_.get(), workspace_.scale_.get(), memspace);
    handler.abs(workspace_.yprev_abs_.get(), workspace_.scale_.get(), memspace);
    handler.max(workspace_.yprev_abs_.get(), workspace_.scale_.get(), workspace_.yprev_abs_.get(), memspace);
    handler.scal(params_.rtol, workspace_.scale_.get(), memspace);
    handler.axpy(1.0, params_.atol.get(), workspace_.scale_.get(), memspace);
    handler.scaleInv(workspace_.scale_.get(), workspace_.out_.get(), memspace);

    return handler.amax(workspace_.out_.get(), memspace);
  }

  template class Rosenbrock<sunrealtype, int>;
} // namespace Integrator