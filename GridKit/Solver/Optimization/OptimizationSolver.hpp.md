```cpp
/**
 * @file OptimizationSolver.hpp
 *
 * @brief Base class for solvers of nonlinear programs.
 *
 * @author Luke Lowery (lukel@tamu.edu)
 */

#pragma once

#include <vector>

#include <GridKit/Solver/Optimization/OptimizationModel.hpp>

namespace GridKit
{
  namespace Optimization
  {
    /**
     * @brief Solves an OptimizationModel, mirroring how the Dynamic
     *        solver family consumes a Model::Evaluator.
     *
     * The solver binds its model at construction and does not own it.
     */
    template <typename scalar_type, typename index_type>
    class OptimizationSolver
    {
    public:
      using ScalarT = scalar_type;
      using IdxT    = index_type;
      using RealT   = typename GridKit::ScalarTraits<ScalarT>::RealT;

      explicit OptimizationSolver(OptimizationModel<ScalarT, IdxT>* model)
        : model_(model)
      {
      }

      virtual ~OptimizationSolver() = default;

      /**
       * @brief Solve the bound model.
       *
       * When @p x is sized to the problem it is the starting point;
       * otherwise the model's default start is used. On success @p x
       * holds the solution.
       *
       * @return 0 on success, positive on a recoverable stop (see the
       *         concrete solver), negative on a hard failure
       */
      [[nodiscard("May fail. Check error code.")]]
      virtual int solve(std::vector<RealT>& x) = 0;

    protected:
      OptimizationModel<ScalarT, IdxT>* model_;
    };
  } // namespace Optimization
} // namespace GridKit
```
