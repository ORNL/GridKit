```cpp
/**
 * @file OptimizationModel.hpp
 *
 * @brief Nonlinear program contract consumed by the Optimization solvers.
 *
 * @author Luke Lowery (lukel@tamu.edu)
 */

#pragma once

#include <GridKit/ScalarTraits.hpp>

namespace GridKit
{
  namespace Optimization
  {
    /**
     * @brief A nonlinear program in the standard form
     *
     *   minimize f(x) subject to c_L <= c(x) <= c_U and x_L <= x <= x_U.
     *
     * Solvers evaluate the problem exclusively through this interface.
     * All arrays are caller-owned, dense, and row-major; the Jacobian is
     * constraintCount() by variableCount(), and the Hessian of the
     * Lagrangian is variableCount() by variableCount(), symmetric, stored
     * full. Unbounded entries use infinity. An equality constraint has
     * equal lower and upper bounds.
     *
     * Every evaluation returns 0 on success, positive when the quantity is
     * unavailable (solvers fall back to an approximation), and negative on
     * a hard failure that aborts the solve.
     */
    template <typename scalar_type, typename index_type>
    class OptimizationModel
    {
    public:
      using ScalarT = scalar_type;
      using IdxT    = index_type;
      using RealT   = typename GridKit::ScalarTraits<ScalarT>::RealT;

      virtual ~OptimizationModel() = default;

      virtual IdxT variableCount() const   = 0;
      virtual IdxT constraintCount() const = 0;

      /**
       * @brief Fill variable bounds, each of length variableCount().
       */
      virtual int variableBounds(RealT* lower, RealT* upper) const = 0;

      /**
       * @brief Fill constraint bounds, each of length constraintCount().
       */
      virtual int constraintBounds(RealT* lower, RealT* upper) const = 0;

      /**
       * @brief Fill the default starting point, length variableCount().
       */
      virtual int startingPoint(RealT* x) const = 0;

      virtual int evaluateObjective(const RealT* x, RealT& value) const = 0;

      /**
       * @brief Fill the objective gradient, length variableCount().
       */
      virtual int evaluateGradient(const RealT* x, RealT* gradient) const = 0;

      /**
       * @brief Fill the constraint values, length constraintCount().
       */
      virtual int evaluateConstraints(const RealT* x, RealT* values) const = 0;

      /**
       * @brief Fill the dense constraint Jacobian, row-major
       *        constraintCount() by variableCount().
       */
      virtual int evaluateJacobian(const RealT* x, RealT* jacobian) const = 0;

      /**
       * @brief Fill the dense Hessian of the Lagrangian,
       *
       *   objective_factor * H(f) + sum over i of multipliers[i] * H(c_i),
       *
       * row-major variableCount() by variableCount(), symmetric. The
       * default declares the Hessian unavailable so solvers use a
       * quasi-Newton approximation.
       */
      virtual int evaluateHessian(const RealT*,
                                  RealT,
                                  const RealT*,
                                  RealT*) const
      {
        return 1;
      }
    };
  } // namespace Optimization
} // namespace GridKit
```
