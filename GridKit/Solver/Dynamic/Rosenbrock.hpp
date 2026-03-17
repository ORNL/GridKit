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

  struct StepControl
  {
    bool   accept;
    double step_size;
  };

  class StepController
  {
  public:
    virtual StepControl nextStep(double err, StepControl prev_step, uint8_t method_order) = 0;
    virtual bool        usesError() const                                                 = 0;
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

      std::string csv_report() const;
      std::string report() const;
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

      std::string report() const;
      Stats&      operator+=(const Stats& other);
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
        return static_cast<bool>(e);
      }

      constexpr bool hasDenseOutput() const
      {
        return static_cast<bool>(H);
      }

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
    [[nodiscard("May fail. Check error code.")]]
    int time_step(double t0, double dt);

    State& error_estimate() const;

    [[nodiscard("May fail. Check error code.")]]
    int calc_dense_coeff();

    [[nodiscard("May fail. Check error code.")]]
    int interp_dense(double theta);
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

    StepControl nextStep(double err, StepControl prev_step, uint8_t method_order) final;

    constexpr bool usesError() const final
    {
      return true;
    }
  };

  class FixedStep : public StepController
  {
    StepControl nextStep(double err, StepControl prev_step, uint8_t method_order) final;

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
