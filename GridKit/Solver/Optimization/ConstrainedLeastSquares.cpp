/**
 * @file ConstrainedLeastSquares.cpp
 *
 * @brief Weighted linear least squares with linear equality constraints.
 *
 * @author Luke Lowery (lukel@tamu.edu)
 */

#include "ConstrainedLeastSquares.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace GridKit
{
  namespace Optimization
  {
    /**
     * @brief Check storage sizes, weight signs, and finiteness.
     */
    template <typename scalar_type, typename index_type>
    int ConstrainedLeastSquares<scalar_type, index_type>::Problem::validate()
        const
    {
      if (variable_count < 1 || residual_count < 1)
      {
        return -1;
      }

      const auto variables  = static_cast<size_t>(variable_count);
      const auto residuals  = static_cast<size_t>(residual_count);
      const auto equalities = static_cast<size_t>(equality_count);

      if (a.size() != residuals * variables || b.size() != residuals)
      {
        return -1;
      }
      if (!weights.empty() && weights.size() != residuals)
      {
        return -1;
      }
      if (a_eq.size() != equalities * variables || b_eq.size() != equalities)
      {
        return -1;
      }
      if (!start.empty() && start.size() != variables)
      {
        return -1;
      }

      for (const auto* values : {&a, &b, &a_eq, &b_eq, &start})
      {
        for (const auto value : *values)
        {
          if (!std::isfinite(value))
          {
            return -2;
          }
        }
      }
      for (const auto value : weights)
      {
        if (!std::isfinite(value) || value < RealT{0})
        {
          return -2;
        }
      }

      return 0;
    }

    /**
     * @brief Bind the problem data and precompute the constant objective
     *        curvature: hessian_ = 2 A^T W A and
     *        gradient_offset_ = -2 A^T W b, so the gradient is
     *        hessian_ x + gradient_offset_.
     */
    template <typename scalar_type, typename index_type>
    ConstrainedLeastSquares<scalar_type, index_type>::ConstrainedLeastSquares(
        Problem problem)
      : problem_(std::move(problem))
    {
      if (problem_.validate() != 0)
      {
        throw std::runtime_error(
            "ConstrainedLeastSquares requires a consistent problem");
      }

      const auto variables = static_cast<size_t>(problem_.variable_count);
      const auto residuals = static_cast<size_t>(problem_.residual_count);

      hessian_.assign(variables * variables, RealT{0});
      gradient_offset_.assign(variables, RealT{0});

      for (size_t m = 0; m < residuals; ++m)
      {
        const RealT  w   = weight(m);
        const RealT* row = &problem_.a[m * variables];

        for (size_t r = 0; r < variables; ++r)
        {
          const RealT scaled = RealT{2} * w * row[r];

          gradient_offset_[r] -= scaled * problem_.b[m];
          for (size_t c = 0; c < variables; ++c)
          {
            hessian_[r * variables + c] += scaled * row[c];
          }
        }
      }
    }

    template <typename scalar_type, typename index_type>
    typename ConstrainedLeastSquares<scalar_type, index_type>::RealT
    ConstrainedLeastSquares<scalar_type, index_type>::weight(
        size_t residual) const
    {
      return problem_.weights.empty() ? RealT{1} : problem_.weights[residual];
    }

    template <typename scalar_type, typename index_type>
    typename ConstrainedLeastSquares<scalar_type, index_type>::IdxT
    ConstrainedLeastSquares<scalar_type, index_type>::variableCount() const
    {
      return problem_.variable_count;
    }

    template <typename scalar_type, typename index_type>
    typename ConstrainedLeastSquares<scalar_type, index_type>::IdxT
    ConstrainedLeastSquares<scalar_type, index_type>::constraintCount() const
    {
      return problem_.equality_count;
    }

    /**
     * @brief Variables are unbounded.
     */
    template <typename scalar_type, typename index_type>
    int ConstrainedLeastSquares<scalar_type, index_type>::variableBounds(
        RealT* lower, RealT* upper) const
    {
      const auto variables = static_cast<size_t>(problem_.variable_count);
      const auto infinity  = std::numeric_limits<RealT>::infinity();

      for (size_t r = 0; r < variables; ++r)
      {
        lower[r] = -infinity;
        upper[r] = infinity;
      }
      return 0;
    }

    /**
     * @brief Equalities: both bounds equal b_eq.
     */
    template <typename scalar_type, typename index_type>
    int ConstrainedLeastSquares<scalar_type, index_type>::constraintBounds(
        RealT* lower, RealT* upper) const
    {
      const auto equalities = static_cast<size_t>(problem_.equality_count);

      for (size_t i = 0; i < equalities; ++i)
      {
        lower[i] = problem_.b_eq[i];
        upper[i] = problem_.b_eq[i];
      }
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int ConstrainedLeastSquares<scalar_type, index_type>::startingPoint(
        RealT* x) const
    {
      const auto variables = static_cast<size_t>(problem_.variable_count);

      for (size_t r = 0; r < variables; ++r)
      {
        x[r] = problem_.start.empty() ? RealT{0} : problem_.start[r];
      }
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int ConstrainedLeastSquares<scalar_type, index_type>::evaluateObjective(
        const RealT* x, RealT& value) const
    {
      const auto variables = static_cast<size_t>(problem_.variable_count);
      const auto residuals = static_cast<size_t>(problem_.residual_count);

      value = RealT{0};
      for (size_t m = 0; m < residuals; ++m)
      {
        const RealT* row = &problem_.a[m * variables];

        RealT residual = -problem_.b[m];
        for (size_t c = 0; c < variables; ++c)
        {
          residual += row[c] * x[c];
        }
        value += weight(m) * residual * residual;
      }
      return 0;
    }

    /**
     * @brief gradient = hessian_ x + gradient_offset_.
     */
    template <typename scalar_type, typename index_type>
    int ConstrainedLeastSquares<scalar_type, index_type>::evaluateGradient(
        const RealT* x, RealT* gradient) const
    {
      const auto variables = static_cast<size_t>(problem_.variable_count);

      for (size_t r = 0; r < variables; ++r)
      {
        RealT value = gradient_offset_[r];
        for (size_t c = 0; c < variables; ++c)
        {
          value += hessian_[r * variables + c] * x[c];
        }
        gradient[r] = value;
      }
      return 0;
    }

    template <typename scalar_type, typename index_type>
    int ConstrainedLeastSquares<scalar_type, index_type>::evaluateConstraints(
        const RealT* x, RealT* values) const
    {
      const auto variables  = static_cast<size_t>(problem_.variable_count);
      const auto equalities = static_cast<size_t>(problem_.equality_count);

      for (size_t i = 0; i < equalities; ++i)
      {
        const RealT* row = &problem_.a_eq[i * variables];

        RealT value = RealT{0};
        for (size_t c = 0; c < variables; ++c)
        {
          value += row[c] * x[c];
        }
        values[i] = value;
      }
      return 0;
    }

    /**
     * @brief The constraint Jacobian is A_eq, constant.
     */
    template <typename scalar_type, typename index_type>
    int ConstrainedLeastSquares<scalar_type, index_type>::evaluateJacobian(
        const RealT*, RealT* jacobian) const
    {
      const auto entries = problem_.a_eq.size();

      for (size_t k = 0; k < entries; ++k)
      {
        jacobian[k] = problem_.a_eq[k];
      }
      return 0;
    }

    /**
     * @brief The Lagrangian Hessian is objective_factor * 2 A^T W A; the
     *        constraints are linear, so the multipliers do not enter.
     */
    template <typename scalar_type, typename index_type>
    int ConstrainedLeastSquares<scalar_type, index_type>::evaluateHessian(
        const RealT*,
        RealT objective_factor,
        const RealT*,
        RealT* hessian) const
    {
      const auto entries = hessian_.size();

      for (size_t k = 0; k < entries; ++k)
      {
        hessian[k] = objective_factor * hessian_[k];
      }
      return 0;
    }

    template class ConstrainedLeastSquares<double, long int>;
    template class ConstrainedLeastSquares<double, size_t>;
  } // namespace Optimization
} // namespace GridKit
