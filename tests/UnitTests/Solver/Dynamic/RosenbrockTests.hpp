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

#include "TrigonometricDae.hpp"
#include <resolve/SystemSolver.hpp>
#include <resolve/vector/Vector.hpp>
#include <resolve/vector/VectorHandler.hpp>
#include <resolve/workspace/LinAlgWorkspaceCpu.hpp>

namespace GridKit
{

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
        model.yp().setToZero(GridKit::memory::HOST);

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
          RealT      error    = 0.0;
          RealT      sol_norm = 0.0;
          const auto solution = model.analyticSolution(final_time);

          // The final solution of the simulation
          const VectorT& state      = model.y();
          const auto*    state_data = state.getData();

          // The difference from the simulated solution to the true solution
          error += std::pow(state_data[0] - solution[0], 2);
          error += std::pow(state_data[1] - solution[1], 2);

          sol_norm += std::pow(solution[0], 2);
          sol_norm += std::pow(solution[1], 2);

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
          model.yp().setToZero(GridKit::memory::HOST);
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
