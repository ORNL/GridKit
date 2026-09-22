#pragma once

#include <cmath>
#include <memory>
#include <stdexcept>

#include <GridKit/LinearAlgebra/Solver/ResolveSystemSolver.hpp>
#include <GridKit/Solver/Dynamic/ImplicitTrapezoidal.hpp>
#include <GridKit/Solver/Dynamic/Native/RmsNorm.hpp>
#include <GridKit/Testing/Testing.hpp>

#include "TrigonometricDae.hpp"
#include <resolve/SystemSolver.hpp>
#include <resolve/workspace/LinAlgWorkspaceCpu.hpp>

namespace GridKit::Testing
{
  class ImplicitTrapezoidalTests
  {
    using Integrator = AnalysisManager::NativeDynamicSolver::ImplicitTrapezoidal<double, int>;
    using Norm       = AnalysisManager::NativeDynamicSolver::RmsNorm<double, int>;
    using State      = LinearAlgebra::Vector<double, int>;

    static Norm::Parameters makeNormParameters()
    {
      auto absolute_tolerance = std::make_unique<State>(Model::TrigonometricDaeEvaluator<double, int>::SIZE);
      if (absolute_tolerance->allocate(memory::HOST) != 0
          || absolute_tolerance->setToConst(1e-12, memory::HOST) != 0)
      {
        throw std::runtime_error("Failed to initialize implicit trapezoidal test tolerances");
      }
      return Norm::Parameters{
          .abs_tol_ = std::move(absolute_tolerance),
          .rel_tol_ = 1e-10,
      };
    }

    struct Fixture
    {
      Fixture()
        : resolve_solver(&workspace, "klu", "klu", "klu"),
          linear_solver(resolve_solver),
          error_norm(makeNormParameters()),
          integrator(&model, linear_solver, vector_handler, &error_norm)
      {
        model.allocate();
        model.initialize();
        resolve_solver.initialize();
        (void) integrator.allocate();
        (void) integrator.initializeSimulation(0.5);
      }

      ReSolve::LinAlgWorkspaceCpu                     workspace;
      Model::TrigonometricDaeEvaluator<double, int>   model;
      ReSolve::SystemSolver                           resolve_solver;
      LinearAlgebra::ResolveSystemSolver<double, int> linear_solver;
      LinearAlgebra::VectorHandler<double, int>       vector_handler;
      Norm                                            error_norm;
      Integrator                                      integrator;
    };

  public:
    TestOutcome differentialAlgebraicSystem()
    {
      TestStatus success        = true;
      double     previous_error = 0.0;

      for (double step_size : {0.05, 0.025, 0.0125})
      {
        Fixture                fixture;
        Integrator::Parameters parameters;
        parameters.step_size_               = step_size;
        parameters.max_steps_               = 200;
        parameters.max_jacobian_refreshes_  = 10;
        parameters.stagnation_ratio_        = 0.1;
        parameters.max_stagnant_iterations_ = 1;

        double y0  = 0.0;
        double y1  = 0.0;
        success   *= fixture.integrator.integrate({2.0}, parameters, [&](double)
                                                {
                                                  y0 = fixture.model.y().getData()[0];
                                                  y1 = fixture.model.y().getData()[1]; })
                   == 0;

        const auto   solution   = fixture.model.analyticSolution(2.0);
        const double error      = std::hypot(y0 - solution[0], y1 - solution[1]);

        if (previous_error != 0.0)
        {
          success *= previous_error / error > 3.9;
        }
        previous_error = error;
      }

      return success.report(__func__);
    }
  };
} // namespace GridKit::Testing
