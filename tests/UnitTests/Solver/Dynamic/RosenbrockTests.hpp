#pragma once

#include <iomanip>
#include <ios>
#include <memory>

#include <GridKit/Model/Evaluator.hpp>
#include <GridKit/Solver/Dynamic/Rosenbrock.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>

#include "GridKit/LinearAlgebra/SparseMatrix/CsrMatrix.hpp"
#include <resolve/vector/Vector.hpp>
#include <resolve/vector/VectorHandler.hpp>
#include <resolve/workspace/LinAlgWorkspaceCpu.hpp>

namespace GridKit
{
  namespace Model
  {
    template <class ScalarT, typename IdxT>
    class TrigonometricDaeEvaluator : public Model::Evaluator<ScalarT, IdxT>
    {
    public:
      using RealT = typename Model::Evaluator<ScalarT, IdxT>::RealT;

      constexpr static size_t SIZE = 2;

      TrigonometricDaeEvaluator()
      {
      }

      int allocate() override
      {
        constexpr size_t NNZ = SIZE * SIZE;

        IdxT*  row_ptrs = new IdxT[SIZE + 1];
        IdxT*  cols     = new IdxT[NNZ];
        RealT* vals     = new RealT[NNZ];

        for (size_t i = 0; i < SIZE + 1; i++)
        {
          row_ptrs[i] = static_cast<IdxT>(i * SIZE);
        }

        for (size_t i = 0; i < SIZE; i++)
        {
          for (size_t j = 0; j < SIZE; j++)
          {
            cols[i * SIZE + j] = static_cast<IdxT>(j);
          }
        }

        csr_jac_ = std::make_unique<GridKit::LinearAlgebra::CsrMatrix<RealT, IdxT>>(SIZE, SIZE, NNZ, &row_ptrs, &cols, &vals);

        return 0;
      }

      int initialize() override
      {
        y_  = {sinh(0.5), tanh(0.5)};
        yp_ = {0.0, 0.0};

        tag_ = {true, false};

        f_ = {0.0, 0.0};
        g_ = {0.0, 0.0};
        return 0;
      }

      IdxT size() override
      {
        return SIZE;
      }

      IdxT nnz() override
      {
        return SIZE * SIZE;
      }

      bool hasJacobian() override
      {
        return true;
      }

      IdxT sizeQuadrature() override
      {
        return 0;
      }

      IdxT sizeParams() override
      {
        return 0;
      }

      int setAbsoluteTolerance([[maybe_unused]] RealT rel_tol) override
      {
        return 0;
      }

      std::vector<ScalarT>& absoluteTolerance() override
      {
        return abs_tol_;
      }

      const std::vector<ScalarT>& absoluteTolerance() const override
      {
        return abs_tol_;
      }

      int tagDifferentiable() override
      {
        return 0;
      }

      int evaluateResidual() override
      {
        double y02 = y_[0] * y_[0];
        double y12 = y_[1] * y_[1];

        f_ = {y02 / (y_[1] * sqrt(std::pow(y_[0] / y_[1], 2) - 1)), //
              y12 + 1 / (1 + y02) - (y02 / y12 - y02)};

        return 0;
      }

      int evaluateJacobian() override
      {
        RealT* vals = csr_jac_->getValues();

        double y1 = y_[0];
        double y2 = y_[1];

        double y12 = y1 * y1;
        double y13 = y12 * y1;
        double y22 = y2 * y2;
        double y23 = y22 * y2;
        double y24 = y22 * y22;

        double tmp  = std::pow(y12 / y22 - 1, 1.5);
        double tmp2 = std::pow(y12 + 1, 2);

        vals[0] = alpha_ + -(-y13 + 2 * y1 * y22) / (y23 * tmp);
        vals[1] = y12 / (y22 * tmp);
        vals[2] = -2 * y1 * (1 / y22 - 1) - (2 * y1) / tmp2;
        vals[3] = (2 * (y12 + y24)) / y23;

        return 0;
      }

      int evaluateIntegrand() override
      {
        return 0;
      }

      int initializeAdjoint() override
      {
        return 0;
      }

      int evaluateAdjointResidual() override
      {
        return 0;
      }

      int evaluateAdjointIntegrand() override
      {
        return 0;
      }

      void updateTime([[maybe_unused]] RealT t, RealT a) override
      {
        alpha_ = a;
      }

      std::vector<ScalarT>& y() override
      {
        return y_;
      }

      const std::vector<ScalarT>& y() const override
      {
        return y_;
      }

      std::vector<ScalarT>& yp() override
      {
        return yp_;
      }

      const std::vector<ScalarT>& yp() const override
      {
        return yp_;
      }

      std::vector<bool>& tag() override
      {
        return tag_;
      }

      const std::vector<bool>& tag() const override
      {
        return tag_;
      }

      std::vector<ScalarT>& yB() override
      {
        return yB_;
      }

      const std::vector<ScalarT>& yB() const override
      {
        return yB_;
      }

      std::vector<ScalarT>& ypB() override
      {
        return ypB_;
      }

      const std::vector<ScalarT>& ypB() const override
      {
        return ypB_;
      }

      std::vector<ScalarT>& param() override
      {
        return param_;
      }

      const std::vector<ScalarT>& param() const override
      {
        return param_;
      }

      std::vector<ScalarT>& param_up() override
      {
        return param_up_;
      }

      const std::vector<ScalarT>& param_up() const override
      {
        return param_up_;
      }

      std::vector<ScalarT>& param_lo() override
      {
        return param_lo_;
      }

      const std::vector<ScalarT>& param_lo() const override
      {
        return param_lo_;
      }

      std::vector<ScalarT>& getResidual() override
      {
        return f_;
      }

      const std::vector<ScalarT>& getResidual() const override
      {
        return f_;
      }

      std::vector<ScalarT>& getIntegrand() override
      {
        return g_;
      }

      const std::vector<ScalarT>& getIntegrand() const override
      {
        return g_;
      }

      std::vector<ScalarT>& getAdjointResidual() override
      {
        return fB_;
      }

      const std::vector<ScalarT>& getAdjointResidual() const override
      {
        return fB_;
      }

      std::vector<ScalarT>& getAdjointIntegrand() override
      {
        return gB_;
      }

      const std::vector<ScalarT>& getAdjointIntegrand() const override
      {
        return gB_;
      }

      IdxT getIDcomponent()
      {
        return 0;
      }

      GridKit::LinearAlgebra::CsrMatrix<RealT, IdxT>* getCsrJacobian() const override
      {
        return csr_jac_.get();
      }

    protected:
      std::vector<ScalarT> y_;
      std::vector<ScalarT> yp_;
      std::vector<bool>    tag_;
      std::vector<ScalarT> f_;
      std::vector<ScalarT> g_;

      std::vector<ScalarT> yB_;
      std::vector<ScalarT> ypB_;
      std::vector<ScalarT> fB_;
      std::vector<ScalarT> gB_;

      std::vector<ScalarT> abs_tol_;

      std::unique_ptr<GridKit::LinearAlgebra::CsrMatrix<RealT, IdxT>> csr_jac_;

      RealT alpha_;

      std::vector<ScalarT> param_;
      std::vector<ScalarT> param_up_;
      std::vector<ScalarT> param_lo_;
    };
  } // namespace Model

  namespace Testing
  {
    template <class ScalarT, typename IdxT>
    class RosenbrockTests
    {
      using Rosenbrock = Integrator::Rosenbrock<ScalarT, IdxT>;

    public:
      TestOutcome test_order(Rosenbrock::Tableau&& tab, double step_exponent_lower, double step_exponent_upper)
      {
        TestStatus success = true;

        uint8_t expected_order = tab.order_;

        Model::TrigonometricDaeEvaluator<ScalarT, IdxT> model;
        model.allocate();
        model.initialize();

        ReSolve::LinAlgWorkspaceCpu                          linear_workspace;
        ReSolve::SystemSolver                                lin_solver(&linear_workspace, "klu", "klu", "klu");
        GridKit::LinearAlgebra::VectorHandler<ScalarT, IdxT> vec_handler;

        lin_solver.initialize();

        Rosenbrock integrator(std::move(tab), &model, lin_solver, vec_handler, nullptr);
        if (integrator.allocate())
        {
          success = false;
          return success.report(__func__);
        }

        Integrator::FixedStep step_controller;

        double              final_time = 2.0;
        std::vector<double> out_times  = {final_time};

        size_t num_samples = 21;

        std::vector<double> step_sizes;
        std::vector<double> errors;

        auto out_cb = [&]([[maybe_unused]] double t)
        {
          double error    = 0.0;
          double sol_norm = 0.0;

          const std::vector<double>& state = model.y();

          error += std::pow(state[0] - sinh(final_time), 2);
          error += std::pow(state[1] - tanh(final_time), 2);

          sol_norm += std::pow(sinh(final_time), 2);
          sol_norm += std::pow(tanh(final_time), 2);

          errors.push_back(std::sqrt(error) / std::sqrt(sol_norm));
        };

        for (size_t i = 0; i < num_samples; i++)
        {
          double step_size = std::pow(10, step_exponent_lower + static_cast<double>(i) * (step_exponent_upper - step_exponent_lower) / static_cast<double>(num_samples - 1));
          double num_steps = round((final_time - 0.5) / step_size);
          step_size        = (final_time - 0.5) / num_steps;
          step_sizes.push_back(step_size);

          model.initialize();
          if (integrator.initializeSimulation(0.5))
          {
            success = false;
            return success.report(__func__);
          }

          typename Rosenbrock::Parameters params;
          params.starting_step_ = step_size;
          params.max_steps_     = static_cast<size_t>(ceil((final_time - 0.5) / step_size)) + 10;
          if (integrator.integrate(out_times, step_controller, params, out_cb))
          {
            success = false;
            return success.report(__func__);
          }
        }

        std::cerr << "Step sizes\n";
        for (double step_size : step_sizes)
        {
          std::cerr << std::scientific << std::setprecision(20) << step_size << ' ';
        }

        std::cerr << "\n\nErrors\n";
        for (double error : errors)
        {
          std::cerr << std::scientific << std::setprecision(20) << error << ' ';
        }
        std::cerr << "\n\n";

        std::vector<double> pairwise_orders;
        double              empirical_order = 0.0;
        for (size_t i = 1; i < num_samples; i++)
        {
          pairwise_orders.push_back((log(errors[i]) - log(errors[i - 1])) / (log(step_sizes[i]) - log(step_sizes[i - 1])));
          empirical_order += pairwise_orders.back();
        }
        empirical_order /= static_cast<double>(num_samples - 1);

        std::cerr << "Empirical order: " << std::fixed << std::setprecision(5) << empirical_order << " v.s. expected " << static_cast<unsigned>(expected_order) << '\n';
        success *= empirical_order > expected_order * 0.85;

        return success.report(__func__);
      }
    };
  } // namespace Testing
} // namespace GridKit
