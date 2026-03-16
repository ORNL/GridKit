#pragma once

#include <cassert>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <ios>
#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

#include <sundials/sundials_nvector.h>

#include "GridKit/Model/Evaluator.hpp"
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

  struct StepControl
  {
    bool   accept;
    double step_size;
  };

  class StepController
  {
  public:
    virtual constexpr StepControl nextStep(double err, StepControl prev_step, uint8_t method_order) = 0;
    virtual constexpr bool        usesError() const                                                 = 0;
  };

  class ErrorNorm
  {
  public:
    virtual double errorNorm(State& err, State& y, State& yprev, ReSolve::VectorHandler& handler, ReSolve::memory::MemorySpace memspace) const = 0;
  };

  template <class ScalarT, typename IdxT>
  class Rosenbrock
  {
    using RealT = typename GridKit::ScalarTraits<ScalarT>::RealT;

  public:
    struct StepInfo
    {
      double sim_time;
      double step_size;
      double next_step_size;
      double err_est;
      size_t step_no;
      bool   skip_lu;
      bool   skip_f;
      bool   accepted;

      std::string csv_report() const
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

      std::string report() const
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
    };

    struct Stats
    {
      std::vector<StepInfo> rejections;
      std::vector<StepInfo> skip_lu_steps;
      size_t                num_steps     = 0;
      size_t                f_evals       = 0;
      size_t                f_skipped     = 0;
      size_t                jac_evals     = 0;
      size_t                decomp_solves = 0;
      double                min_step      = INFINITY;
      double                max_step      = 0;

      std::string report() const
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

      Stats& operator+=(const Stats& other)
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
    };

    struct Parameters
    {
      double atol          = 1e-5;
      double rtol          = 1e-5;
      double starting_step = 1e-5;
      size_t max_steps     = 2000;
      bool   skip_lu       = false;
    };

    struct Tableau
    {

      size_t num_stages;

      /// @brief The coefficient along the diagonal of the Gamma matrix.
      RealT gamma;

      /// @brief A vector of sums of rows of the alpha matrix. These are the classic
      ///        Runge-Kutta 'c' coefficients, or abscissae.
      std::unique_ptr<RealT[]> alpha_sum;
      /// @brief A vector of sums of rows of the Gamma matrix.
      std::unique_ptr<RealT[]> gamma_sum;
      /// @brief A vector of weights for constructing the final solution from the stages
      std::unique_ptr<RealT[]> m;

      /// @brief Coefficients for the embedded method. If `HAS_EMBEDDED == false`, then
      ///        this is nothing.
      std::unique_ptr<RealT[]> e;

      /// @brief The transformed A coefficient matrix. Strictly lower triangular.
      std::unique_ptr<RealT[]> A;
      /// @brief The transformed C coefficient matrix. Strictly lower triangular.
      std::unique_ptr<RealT[]> C;

      /// @brief The matrix of dense coefficients.
      std::unique_ptr<RealT[]> H;

      /// @brief What ODE order these coefficients satisfy.
      uint8_t order;

      /// @brief Whether or not these coefficients satisfy Rosenbrock-Krylov (ROK) order conditions up to
      /// `order`.
      bool is_krylov;
      /// @brief Whether or not these coefficients satisfy Rosenbrock-W (ROW) order conditions up to `order`.
      bool is_w;
      /// @brief Whether or not these coefficients satisfy DAE order conditions up to `order`.
      bool is_dae;

      /// @brief Whether or not this tableau contains an embedded method
      constexpr bool has_embedded() const
      {
        return e == nullptr;
      }

      constexpr bool hasDenseOutput() const
      {
        return static_cast<bool>(H);
      }

      constexpr RealT getA(size_t row, size_t col) const
      {
        return A[row * num_stages + col];
      }

      constexpr bool can_reuse_asum(size_t stage) const
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

      constexpr bool can_reuse_asum_for_out() const
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

      constexpr std::tuple<bool, size_t> error_estimator_stage() const
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

      static constexpr Tableau lin_implicit_euler()
      {
        constexpr size_t num_stages = 1;

        Tableau re = {
            .num_stages = 1,
            .gamma      = 1.0,
            .alpha_sum  = std::make_unique_for_overwrite<RealT[]>(num_stages),
            .gamma_sum  = std::make_unique_for_overwrite<RealT[]>(num_stages),
            .m          = std::make_unique_for_overwrite<RealT[]>(num_stages),
            .e          = {},
            .A          = std::make_unique_for_overwrite<RealT[]>(num_stages * num_stages),
            .C          = std::make_unique_for_overwrite<RealT[]>(num_stages * num_stages),
            .H          = {},
            .order      = 1,
            .is_krylov  = true,
            .is_w       = false,
            .is_dae     = true,
        };

        re.alpha_sum[0] = 0.0;

        re.gamma_sum[0] = 1.0;

        re.m[0] = 1.0;

        re.A[0] = 0.0;

        re.C[0] = 2.0;

        return re;
      }
    };

  private:
    double step_size_      = 0;
    double prev_step_size_ = 0;

    bool skip_lu_ = false;
    bool skip_f_  = false;

    /// @brief Keeps track of whether or not the integrator currently has valid dense coefficients.
    ///        i.e. they have been computed and haven't been invalidated by taking another step.
    bool dense_coefficients_valid = false;

    Tableau tab_;

    GridKit::Model::Evaluator<ScalarT, IdxT>* model_;
    ReSolve::SystemSolver&                    lin_solver_;
    ReSolve::VectorHandler&                   vector_handler_;
    const ErrorNorm*                          err_norm_;
    ReSolve::memory::MemorySpace              memspace_;

    double                 current_time_;
    std::unique_ptr<State> y_prev_;
    std::unique_ptr<State> y_cur_;
    std::unique_ptr<State> y_interp_;

    Parameters params_;
    Stats      stats_;

    Rosenbrock(const Rosenbrock& other) = delete;
    Rosenbrock(Rosenbrock&& other)      = delete;

  public:
    Rosenbrock(Tableau&&                                 tab,
               GridKit::Model::Evaluator<ScalarT, IdxT>* model,
               ReSolve::SystemSolver&                    lin_solver,
               ReSolve::VectorHandler&                   vector_handler,
               const ErrorNorm*                          err_norm,
               ReSolve::memory::MemorySpace              memspace = ReSolve::memory::HOST)
      : tab_(std::move(tab)),
        model_(model),
        lin_solver_(lin_solver),
        vector_handler_(vector_handler),
        err_norm_(err_norm),
        memspace_(memspace)
    {
    }

    int allocate()
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

      y_prev_->allocate(memspace_);
      y_cur_->allocate(memspace_);
      y_interp_->allocate(memspace_);
      asum_->allocate(memspace_);
      csum_->allocate(memspace_);
      RHS_->allocate(memspace_);
      RHS_first_stage_->allocate(memspace_);
      dFdt_->allocate(memspace_);
      mass_->allocate(memspace_);
      y_new_->allocate(memspace_);

      if (tab_.e)
      {
        auto [can_use_stage, err_stage] = tab_.error_estimator_stage();
        if (!can_use_stage)
        {
          err_est_ = std::make_unique<State>(size);
          err_est_->allocate(memspace_);
        }
      }

      stages_ = std::make_unique<std::unique_ptr<State>[]>(tab_.num_stages);
      for (size_t i = 0; i < tab_.num_stages; i++)
      {
        stages_[i] = std::make_unique<State>(size);
        stages_[i]->allocate(memspace_);
      }

      if (tab_.order > 2)
      {
        dense_coeff_ = std::make_unique<std::unique_ptr<State>[]>(tab_.order - 2);
        for (size_t i = 0; i < tab_.order - 2; i++)
        {
          dense_coeff_[i] = std::make_unique<State>(size);
          dense_coeff_[i]->allocate(memspace_);
        }
      }

      jacobian_ = std::make_unique<ReSolve::matrix::Csr>(size, size, model_->getCsrJacobian()->getNnz());

      return 0;
    }

    int initializeSimulation(RealT t0)
    {
      current_time_ = t0;
      y_cur_->copyFromExternal(model_->y().data(), memspace_, memspace_);
      jacobian_analyzed_ = false;

      GridKit::LinearAlgebra::CsrMatrix<RealT, IdxT>* model_jacobian = model_->getCsrJacobian();
      jacobian_->setDataPointers(
          model_jacobian->getRowData(),
          model_jacobian->getColData(),
          model_jacobian->getValues(),
          memspace_);
      lin_solver_.setMatrix(jacobian_.get());
      lin_solver_.analyze();
      lin_solver_.preconditionerSetup();

      assert(model_->tag().size() == static_cast<size_t>(model_->size()));
      std::unique_ptr<RealT[]> mass = std::make_unique<RealT[]>(model_->tag().size());
      for (size_t i = 0; i < static_cast<size_t>(model_->size()); i++)
      {
        mass[i] = model_->tag()[i] ? 1.0 : 0.0;
      }
      mass_->copyFromExternal(mass.get(), memspace_, memspace_);

      stats_ = Stats();

      return 0;
    }

    int integrate(const std::vector<double>&                          out_times,
                  StepController&                                     step_controller,
                  Parameters                                          params  = {},
                  std::optional<std::function<void(double)>>          out_cb  = {},
                  std::optional<std::function<void(const StepInfo&)>> step_cb = {})
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
          // TODO: This assumes the step cannot fail.
          time_step(current_time_, step_size_);

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
            y_cur_->copyToExternal(model_->y().data(), memspace_, memspace_);
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
              calc_dense_coeff();
              dense_coefficients_valid = true;
            }

            double theta = (out_time - current_time_) / prev_step_size_ + 1;
            interp_dense(theta);
          }
          else
          {
            // TODO: Put code for alternative interpolation (Abdou) here
            double theta = (out_time - current_time_) / prev_step_size_ + 1;
            y_interp_->copyFromExternal(y_prev_.get(), memspace_, memspace_);
            vector_handler_.scal(1 - theta, y_interp_.get(), memspace_);
            vector_handler_.axpy(theta, y_cur_.get(), y_interp_.get(), memspace_);
          }

          if (out_cb)
          {
            y_interp_->copyToExternal(model_->y().data(), memspace_, memspace_);
            model_->updateTime(out_time, 0.0);

            (*out_cb)(out_time);
          }
        }
        else
        {
          y_interp_->copyFromExternal(y_cur_.get(), memspace_, memspace_);
          break;
        }
      }

      return 0;
    }

  private:
    std::unique_ptr<State> asum_;
    std::unique_ptr<State> csum_;
    std::unique_ptr<State> RHS_;
    std::unique_ptr<State> RHS_first_stage_;
    std::unique_ptr<State> dFdt_;
    std::unique_ptr<State> mass_;
    std::unique_ptr<State> y_new_;
    std::unique_ptr<State> err_est_;

    std::unique_ptr<ReSolve::matrix::Csr> jacobian_;
    bool                                  jacobian_analyzed_ = false;

    std::unique_ptr<std::unique_ptr<State>[]> stages_;
    std::unique_ptr<std::unique_ptr<State>[]> dense_coeff_;

  public:
    int time_step(double t0, double dt)
    {
      bool y0_copied = false;

      // Form the left-hand side of the system
      // Can sometimes be skipped if the method is a w-method
      [[likely]]
      if (!tab_.is_w || !skip_lu_)
      {
        y_cur_->copyToExternal(model_->y().data(), memspace_, memspace_);
        y0_copied = true;
        model_->updateTime(t0, -1.0 / (dt * tab_.gamma));
        model_->evaluateJacobian();
        GridKit::LinearAlgebra::CsrMatrix<RealT, IdxT>* model_jacobian = model_->getCsrJacobian();
        jacobian_->setDataPointers(
            model_jacobian->getRowData(),
            model_jacobian->getColData(),
            model_jacobian->getValues(),
            memspace_);

        [[likely]]
        if (jacobian_analyzed_)
        {
          lin_solver_.refactorize();
        }
        else
        {
          lin_solver_.factorize();
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
          y_cur_->copyToExternal(model_->y().data(), memspace_, memspace_);
          y0_copied = true;
        }
        model_->updateTime(t0, 0.0);
        model_->evaluateResidual();
        RHS_first_stage_->copyFromExternal(model_->getResidual().data(), memspace_, memspace_);
        vector_handler_.scal(-1, RHS_first_stage_.get(), memspace_);

        stats_.f_evals++;
      }
      lin_solver_.solve(RHS_first_stage_.get(), stages_[0].get());
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
          asum_->copyFromExternal(y_cur_.get(), memspace_, memspace_);

          for (size_t j = 0; j < i; j++)
          {
            if (tab_.A[tab_.num_stages * i + j] != 0.0)
              vector_handler_.axpy(tab_.A[tab_.num_stages * i + j], stages_[j].get(), asum_.get(), memspace_);
          }
        }

        // Calculate csum
        // TODO: Since csum is multiplied by the mass matrix, we can reduce calculations by just not calculating some indices
        csum_->setToZero(memspace_);
        for (size_t j = 0; j < i; j++)
        {
          if (tab_.C[i * tab_.num_stages + j] != 0.0)
          {
            vector_handler_.axpy(tab_.C[i * tab_.num_stages + j] / dt, stages_[j].get(), csum_.get(), memspace_);
          }
        }

        // TODO: non-autonomous model
        asum_->copyToExternal(model_->y().data(), memspace_, memspace_);
        model_->updateTime(t0 + tab_.alpha_sum[i] * dt, 0.0);
        model_->evaluateResidual();
        RHS_->copyFromExternal(model_->getResidual().data(), memspace_, memspace_);

        vector_handler_.scal(-1, RHS_.get(), memspace_);
        vector_handler_.scal(mass_.get(), csum_.get(), memspace_);
        vector_handler_.axpy(-1, csum_.get(), RHS_.get(), memspace_);

        lin_solver_.solve(RHS_.get(), stages_[i].get());
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
        y_new_->copyFromExternal(y_cur_.get(), memspace_, memspace_);

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

    State& error_estimate() const
    {
      auto [can_use_stage, err_stage] = tab_.error_estimator_stage();

      if (can_use_stage)
      {
        return *stages_[err_stage];
      }
      else
      {
        err_est_->copyFromExternal(stages_[0].get(), memspace_, memspace_);
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

    constexpr bool has_dense_output()
    {
      return static_cast<bool>(tab_.H);
    }

    void calc_dense_coeff()
    {
      if (tab_.order > 2)
      {
        for (size_t j = 0; j < tab_.order - 2; j++)
        {
          dense_coeff_[j]->setToZero(memspace_);
        }

        for (size_t i = 0; i < tab_.num_stages; i++)
        {
          for (size_t j = 0; j < tab_.order - 2; j++)
          {
            vector_handler_.axpy(tab_.H[j * (tab_.num_stages - 2) + i], stages_[i].get(), dense_coeff_[j].get(), memspace_);
          }
        }
      }
    }

    // TODO: Maybe this can be integrated into OneStepIntegrator?
    int interp_dense(double theta)
    {
      double one = 1.0;
      if (tab_.order > 2)
      {
        y_interp_->copyFromExternal(dense_coeff_[tab_.order - 3].get(), memspace_, memspace_);

        for (size_t i = 1; i < tab_.order - 2; i++)
        {
          vector_handler_.scal(theta, y_interp_.get(), memspace_);
          vector_handler_.axpy(one, dense_coeff_[tab_.order - 3 - i].get(), y_interp_.get(), memspace_);
        }
      }
      else
      {
        y_interp_->setToZero(memspace_);
      }

      double omt = 1 - theta;

      vector_handler_.scal(omt, y_interp_.get(), memspace_);
      vector_handler_.axpy(one, y_cur_.get(), y_interp_.get(), memspace_);
      vector_handler_.scal(theta, y_interp_.get(), memspace_);
      vector_handler_.axpy(omt, y_prev_.get(), y_interp_.get(), memspace_);

      return 0;
    }
  };

  class AdaptiveStep : public StepController
  {
    struct Parameters
    {
      double facmin   = 0.2;
      double facmax   = 5.0;
      double facscale = 0.9;
    } params_;

  public:
    AdaptiveStep(const Parameters& params) : params_(params)
    {
    }

    constexpr StepControl nextStep(double err, StepControl prev_step, uint8_t method_order) final
    {
      StepControl next_step = prev_step;

      double h_mult = std::min(params_.facmax, std::max(params_.facscale * std::pow(err, -1.0 / method_order), params_.facmin));

      next_step.accept     = err <= 1;
      next_step.step_size *= h_mult;

      return next_step;
    }

    constexpr bool usesError() const final
    {
      return true;
    }
  };

  class FixedStep : public StepController
  {
    StepControl nextStep([[maybe_unused]] double err, [[maybe_unused]] StepControl prev_step, [[maybe_unused]] uint8_t method_order)
    {
      return StepControl{
          .accept    = true,
          .step_size = prev_step.step_size,
      };
    }

    bool usesError() const final
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

    double errorNorm(State& err, State& y, State& yprev, ReSolve::VectorHandler& handler, ReSolve::memory::MemorySpace memspace) const final
    {
      double one = 1.0;
      workspace_.out_->copyFromExternal(&err, memspace, memspace);
      workspace_.scale_->copyFromExternal(&y, memspace, memspace);
      workspace_.yprev_abs_->copyFromExternal(&yprev, memspace, memspace);

      handler.abs(workspace_.scale_.get(), workspace_.scale_.get(), memspace);
      handler.abs(workspace_.yprev_abs_.get(), workspace_.scale_.get(), memspace);
      handler.max(workspace_.yprev_abs_.get(), workspace_.scale_.get(), workspace_.yprev_abs_.get(), memspace);
      handler.scal(params_.rtol, workspace_.scale_.get(), memspace);
      handler.axpy(one, params_.atol.get(), workspace_.scale_.get(), memspace);
      handler.scaleInv(workspace_.scale_.get(), workspace_.out_.get(), memspace);

      return handler.amax(workspace_.out_.get(), memspace);
    }
  };
} // namespace Integrator