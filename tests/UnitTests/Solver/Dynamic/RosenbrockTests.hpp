#pragma once

#include <iomanip>
#include <ios>
#include <memory>

#include <GridKit/LinearAlgebra/Solver/ResolveSystemSolver.hpp>
#include <GridKit/LinearAlgebra/SparseMatrix/CsrMatrix.hpp>
#include <GridKit/Model/Evaluator.hpp>
#include <GridKit/Solver/Dynamic/Native/FixedStep.hpp>
#include <GridKit/Solver/Dynamic/Rosenbrock.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>

#include <resolve/SystemSolver.hpp>
#include <resolve/vector/Vector.hpp>
#include <resolve/vector/VectorHandler.hpp>
#include <resolve/workspace/LinAlgWorkspaceCpu.hpp>

namespace GridKit
{
  namespace Model
  {
    /**
     * @brief A test DAE (referred to as the "trigonometric DAE") that is useful for evaluating order
     * of DAE integrators. The system is 2-dimensional, with one differential variable and one algebraic variable.
     * The DAE is index 1 and the derivatives of the model are non-vanishing.
     *
     * The valid simulation time interval is \f([0.5,2]\f) and the model is initialized at \f(t = 0.5\f).
     *
     */
    template <class ScalarT, typename IdxT>
    class TrigonometricDaeEvaluator : public Model::Evaluator<ScalarT, IdxT>
    {
    public:
      using RealT   = typename Model::Evaluator<ScalarT, IdxT>::RealT;
      using VectorT = typename Model::Evaluator<ScalarT, IdxT>::VectorT;

      constexpr static size_t SIZE = 2;

      TrigonometricDaeEvaluator()
      {
      }

      int allocate() override
      {
        constexpr size_t NNZ = SIZE * SIZE;

        y_.resize(SIZE);
        yp_.resize(SIZE);
        f_.resize(SIZE);

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
        auto* y  = y_.getData();
        auto* yp = yp_.getData();
        auto* f  = f_.getData();

        y[0]  = sinh(0.5);
        y[1]  = tanh(0.5);
        yp[0] = 0.0;
        yp[1] = 0.0;

        tag_ = {true, false};

        f[0] = 0.0;
        f[1] = 0.0;

        y_.setDataUpdated();
        yp_.setDataUpdated();
        f_.setDataUpdated();

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

      VectorT& absoluteTolerance() override
      {
        return abs_tol_;
      }

      const VectorT& absoluteTolerance() const override
      {
        return abs_tol_;
      }

      int tagDifferentiable() override
      {
        return 0;
      }

      int evaluateResidual() override
      {
        const auto* y = y_.getData();
        auto*       f = f_.getData();

        ScalarT y02 = y[0] * y[0];
        ScalarT y12 = y[1] * y[1];

        f[0] = y02 / (y[1] * std::sqrt(std::pow(y[0] / y[1], 2) - 1));
        f[1] = y12 + 1 / (1 + y02) - (y02 / y12 - y02);

        f_.setDataUpdated();

        return 0;
      }

      int evaluateJacobian() override
      {
        RealT* vals = csr_jac_->getValues();

        const auto* y = y_.getData();

        ScalarT y1 = y[0];
        ScalarT y2 = y[1];

        ScalarT y12 = y1 * y1;
        ScalarT y13 = y12 * y1;
        ScalarT y22 = y2 * y2;
        ScalarT y23 = y22 * y2;
        ScalarT y24 = y22 * y22;

        ScalarT tmp  = std::pow(y12 / y22 - 1, 1.5);
        ScalarT tmp2 = std::pow(y12 + 1, 2);

        vals[0] = static_cast<RealT>(alpha_ + -(-y13 + 2 * y1 * y22) / (y23 * tmp));
        vals[1] = static_cast<RealT>(y12 / (y22 * tmp));
        vals[2] = static_cast<RealT>(-2 * y1 * (1 / y22 - 1) - (2 * y1) / tmp2);
        vals[3] = static_cast<RealT>((2 * (y12 + y24)) / y23);

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

      VectorT& y() override
      {
        return y_;
      }

      const VectorT& y() const override
      {
        return y_;
      }

      VectorT& yp() override
      {
        return yp_;
      }

      const VectorT& yp() const override
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

      VectorT& yB() override
      {
        return yB_;
      }

      const VectorT& yB() const override
      {
        return yB_;
      }

      VectorT& ypB() override
      {
        return ypB_;
      }

      const VectorT& ypB() const override
      {
        return ypB_;
      }

      VectorT& param() override
      {
        return param_;
      }

      const VectorT& param() const override
      {
        return param_;
      }

      VectorT& param_up() override
      {
        return param_up_;
      }

      const VectorT& param_up() const override
      {
        return param_up_;
      }

      VectorT& param_lo() override
      {
        return param_lo_;
      }

      const VectorT& param_lo() const override
      {
        return param_lo_;
      }

      VectorT& getResidual() override
      {
        return f_;
      }

      const VectorT& getResidual() const override
      {
        return f_;
      }

      VectorT& getIntegrand() override
      {
        return g_;
      }

      const VectorT& getIntegrand() const override
      {
        return g_;
      }

      VectorT& getAdjointResidual() override
      {
        return fB_;
      }

      const VectorT& getAdjointResidual() const override
      {
        return fB_;
      }

      VectorT& getAdjointIntegrand() override
      {
        return gB_;
      }

      const VectorT& getAdjointIntegrand() const override
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
      VectorT           y_;
      VectorT           yp_;
      std::vector<bool> tag_;
      VectorT           f_;
      VectorT           g_;

      VectorT yB_;
      VectorT ypB_;
      VectorT fB_;
      VectorT gB_;

      VectorT abs_tol_;

      std::unique_ptr<GridKit::LinearAlgebra::CsrMatrix<RealT, IdxT>> csr_jac_;

      RealT alpha_;

      VectorT param_;
      VectorT param_up_;
      VectorT param_lo_;
    };
  } // namespace Model

  namespace Testing
  {
    template <class ScalarT, typename IdxT>
    class RosenbrockTests
    {
      using Rosenbrock = AnalysisManager::NativeDynamicSolver::Rosenbrock<ScalarT, IdxT>;
      using RealT      = typename GridKit::ScalarTraits<ScalarT>::RealT;
      using VectorT    = typename Model::Evaluator<ScalarT, IdxT>::VectorT;

    public:
      /**
       * @brief Test a Rosenbrock tableau by verifying its order empirically. Empirical order is calculated
       * by running a convergence test on the problem modeled in \ref Model::TrigonometricDaeEvaluator. 21
       * simulations are run on the problem with a fixed-step step size controller with a number of steps
       * logarithmically distributed in \f(\left[10^a, 10^b\right]\f) where \f(a\f) is `step_exponent_lower`
       * and \f(b\f) is `step_exponent_upper`. The average slope between pairs of simulations is taken to be
       * the empirical order. The test succeeds if the empirical order is at least 85% of the expected theoretical order.
       *
       * @param tab The tableau to test
       * @param step_exponent_lower The exponent describing the smallest number of steps to take during a simulation.
       * @param step_exponent_upper The exponent describing the largest number of steps to take during a simulation.
       */
      TestOutcome test_order(Rosenbrock::Tableau&& tab, RealT step_exponent_lower, RealT step_exponent_upper)
      {
        TestStatus success = true;

        // Tableaus keep track of their theoretical order. We will attempt to match the empirical order to this.
        uint8_t expected_order = tab.order_;

        // Setup the model, linear solver, and integrator
        Model::TrigonometricDaeEvaluator<ScalarT, IdxT> model;
        model.allocate();
        model.initialize();

        ReSolve::LinAlgWorkspaceCpu                                linear_workspace;
        ReSolve::SystemSolver                                      resolve_solver(&linear_workspace, "klu", "klu", "klu");
        GridKit::LinearAlgebra::VectorHandler<ScalarT, IdxT>       vec_handler;
        GridKit::LinearAlgebra::ResolveSystemSolver<ScalarT, IdxT> lin_solver(resolve_solver);

        resolve_solver.initialize();

        Rosenbrock integrator(std::move(tab), &model, lin_solver, vec_handler, nullptr);
        if (integrator.allocate())
        {
          success = false;
          return success.report(__func__);
        }

        AnalysisManager::NativeDynamicSolver::FixedStep<RealT> step_controller;

        // The number of simulations to run to calculate empirical order. Must be at least 2,
        // since empirical order is calculated pairwise between simulations.
        size_t num_samples = 21;

        // Output vectors to keep track of data used to calculate empirical order.
        // Each simulation will record its step size and final error.
        std::vector<RealT> step_sizes;
        std::vector<RealT> errors;

        // An output callback to populate step_sizes and errors. Each simulation will call this callback
        // at the end of the simulation (final_time). It will then calculate the final error by comparing the
        // solution of the simulation to the true answer.
        RealT              final_time = 2.0;
        std::vector<RealT> out_times  = {final_time};
        auto               out_cb     = [&]([[maybe_unused]] RealT t)
        {
          RealT error    = 0.0;
          RealT sol_norm = 0.0;

          // The final solution of the simulation
          const VectorT& state      = model.y();
          const auto*    state_data = state.getData();

          // The difference from the simulated solution to the true solution
          error += std::pow(state_data[0] - sinh(final_time), 2);
          error += std::pow(state_data[1] - tanh(final_time), 2);

          sol_norm += std::pow(sinh(final_time), 2);
          sol_norm += std::pow(tanh(final_time), 2);

          // Error relative to the true solution
          errors.push_back(std::sqrt(error) / std::sqrt(sol_norm));
        };

        // Perform all of the simulations and populate step_sizes and errors using out_cb
        for (size_t i = 0; i < num_samples; i++)
        {
          // Logarithmically distribute step_size based on step_exponent_lower and step_exponent_upper.
          // Round it to ensure num_steps is integral (do not invoke the dense output, which can add additional errors).
          RealT step_size = std::pow(10, step_exponent_lower + static_cast<RealT>(i) * (step_exponent_upper - step_exponent_lower) / static_cast<RealT>(num_samples - 1));
          RealT num_steps = round((final_time - 0.5) / step_size);
          step_size       = (final_time - 0.5) / num_steps;
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

        // Print output data
        std::cout << "Step sizes\n";
        for (RealT step_size : step_sizes)
        {
          std::cout << std::scientific << std::setprecision(20) << step_size << "\n";
        }

        std::cout << "\nErrors\n";
        for (RealT error : errors)
        {
          std::cout << std::scientific << std::setprecision(20) << error << "\n";
        }
        std::cout << "\n";

        // Calculate empirical order. Each pairwise order is calculated, then averaged.
        std::vector<RealT> pairwise_orders;
        RealT              empirical_order = 0.0;
        for (size_t i = 1; i < num_samples; i++)
        {
          pairwise_orders.push_back((log(errors[i]) - log(errors[i - 1])) / (log(step_sizes[i]) - log(step_sizes[i - 1])));
          empirical_order += pairwise_orders.back();
        }
        empirical_order /= static_cast<RealT>(num_samples - 1);

        // Print test result - observed empirical order and the expected order.
        std::cout << "Empirical order: " << std::fixed << std::setprecision(5) << empirical_order << "\n"
                  << "Expected order: " << static_cast<unsigned>(expected_order) << "\n";
        success *= empirical_order > expected_order * 0.85;

        return success.report(__func__);
      }
    };
  } // namespace Testing
} // namespace GridKit
