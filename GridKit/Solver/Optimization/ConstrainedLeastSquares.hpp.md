```cpp
/**
 * @file ConstrainedLeastSquares.hpp
 *
 * @brief Weighted linear least squares with linear equality constraints.
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
     * @brief The convex program
     *
     *   minimize sum over m of weights[m] * (A x - b)_m^2
     *   subject to A_eq x = b_eq,
     *
     * posed as an OptimizationModel. The objective Hessian 2 A^T W A and
     * the gradient offset -2 A^T W b are precomputed at construction, so
     * the exact Hessian is available to the solver at no per-iteration
     * cost. Both estimator quadratic programs of vector fitting are
     * instances of this model.
     */
    template <typename scalar_type, typename index_type>
    class ConstrainedLeastSquares final
      : public OptimizationModel<scalar_type, index_type>
    {
    public:
      using ScalarT = scalar_type;
      using IdxT    = index_type;
      using RealT   = typename GridKit::ScalarTraits<ScalarT>::RealT;

      /// Dense problem data, row-major.
      struct Problem
      {
        std::vector<RealT> a;       ///< residual_count x variable_count.
        std::vector<RealT> b;       ///< residual_count.
        std::vector<RealT> weights; ///< residual_count, empty is uniform.
        std::vector<RealT> a_eq;    ///< equality_count x variable_count.
        std::vector<RealT> b_eq;    ///< equality_count.
        std::vector<RealT> start;   ///< variable_count, empty is zeros.
        IdxT               variable_count{0};
        IdxT               residual_count{0};
        IdxT               equality_count{0};

        /**
         * @brief Check storage sizes, weight signs, and finiteness.
         *
         * @return 0 when valid, -1 on inconsistent sizes, -2 on a
         *         nonfinite value or negative weight
         */
        [[nodiscard("May fail. Check error code.")]]
        int validate() const;
      };

      /**
       * @brief Bind the problem data.
       *
       * @pre problem.validate() returns 0; the constructor throws
       *      std::runtime_error otherwise
       */
      explicit ConstrainedLeastSquares(Problem problem);

      IdxT variableCount() const override;
      IdxT constraintCount() const override;

      int variableBounds(RealT* lower, RealT* upper) const override;
      int constraintBounds(RealT* lower, RealT* upper) const override;
      int startingPoint(RealT* x) const override;

      int evaluateObjective(const RealT* x, RealT& value) const override;
      int evaluateGradient(const RealT* x, RealT* gradient) const override;
      int evaluateConstraints(const RealT* x, RealT* values) const override;
      int evaluateJacobian(const RealT* x, RealT* jacobian) const override;
      int evaluateHessian(const RealT* x,
                          RealT        objective_factor,
                          const RealT* multipliers,
                          RealT*       hessian) const override;

    private:
      RealT weight(size_t residual) const;

      Problem            problem_;
      std::vector<RealT> hessian_;         ///< 2 A^T W A, row-major.
      std::vector<RealT> gradient_offset_; ///< -2 A^T W b.
    };
  } // namespace Optimization
} // namespace GridKit
```
