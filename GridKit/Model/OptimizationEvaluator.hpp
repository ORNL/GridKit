/**
 * @file OptimizationEvaluator.hpp
 * @brief Abstract nonlinear program with exact sparse derivatives.
 */

#pragma once

#include <memory>
#include <vector>

#include <GridKit/LinearAlgebra/SparseMatrix/CsrMatrix.hpp>
#include <GridKit/LinearAlgebra/Vector/Vector.hpp>
#include <GridKit/ScalarTraits.hpp>

namespace GridKit
{
  namespace Model
  {
    /**
     * @brief Minimize \f$f(x)\f$ subject to \f$g_l \le g(x) \le g_u\f$ and
     * \f$x_l \le x \le x_u\f$
     *
     * Evaluations are at `x()` and `lambda()`, and return 0 on success.
     */
    template <typename scalar_type, typename index_type>
    class OptimizationEvaluator
    {
    public:
      using ScalarT     = scalar_type;
      using IdxT        = index_type;
      using RealT       = typename ScalarTraits<ScalarT>::RealT;
      using VectorT     = LinearAlgebra::Vector<ScalarT, IdxT>;
      using RealVectorT = LinearAlgebra::Vector<RealT, IdxT>;
      using CsrMatrixT  = LinearAlgebra::CsrMatrix<RealT, IdxT>;

      virtual ~OptimizationEvaluator() = default;

      /// Size the problem, tag the inequalities, and set the sparsity patterns
      virtual int allocate() = 0;

      /// Fill `tag()`, true for an inequality and false for an equality
      virtual int tagInequality() = 0;

      /// Set the bounds and the flat start
      virtual int initialize() = 0;

      /// Evaluate \f$f\f$ into `f()`
      virtual int evaluateObjective() = 0;

      /// Evaluate \f$\nabla f\f$ into `gradient()`
      virtual int evaluateGradient() = 0;

      /// Evaluate \f$g\f$ into `g()`
      virtual int evaluateConstraints() = 0;

      /// Evaluate the Jacobian of \f$g\f$ into `getCsrJacobian()`
      virtual int evaluateJacobian() = 0;

      /// Evaluate the lower Hessian of \f$\sigma f + \lambda^T g\f$ into `getCsrHessian()`
      virtual int evaluateHessian(RealT sigma) = 0;

      /// Number of variables
      IdxT size() const
      {
        return x_.getSize();
      }

      /// Number of constraints
      IdxT sizeConstraints() const
      {
        return g_.getSize();
      }

      /// Objective
      const ScalarT& f() const
      {
        return f_;
      }

      /// Variables
      VectorT& x()
      {
        return x_;
      }

      const VectorT& x() const
      {
        return x_;
      }

      /// Variable lower bounds
      const RealVectorT& xLower() const
      {
        return x_lower_;
      }

      /// Variable upper bounds
      const RealVectorT& xUpper() const
      {
        return x_upper_;
      }

      /// Objective gradient
      const RealVectorT& gradient() const
      {
        return gradient_;
      }

      /// Constraints
      const VectorT& g() const
      {
        return g_;
      }

      /// Constraint lower bounds
      const RealVectorT& gLower() const
      {
        return g_lower_;
      }

      /// Constraint upper bounds
      const RealVectorT& gUpper() const
      {
        return g_upper_;
      }

      /// True for each inequality
      const std::vector<bool>& tag() const
      {
        return tag_;
      }

      /// Constraint multipliers
      RealVectorT& lambda()
      {
        return lambda_;
      }

      const RealVectorT& lambda() const
      {
        return lambda_;
      }

      /// Lower-bound multipliers
      RealVectorT& zLower()
      {
        return z_lower_;
      }

      const RealVectorT& zLower() const
      {
        return z_lower_;
      }

      /// Upper-bound multipliers
      RealVectorT& zUpper()
      {
        return z_upper_;
      }

      const RealVectorT& zUpper() const
      {
        return z_upper_;
      }

      /// Constraint Jacobian
      CsrMatrixT* getCsrJacobian() const
      {
        return jacobian_.get();
      }

      /// Lower Lagrangian Hessian
      CsrMatrixT* getCsrHessian() const
      {
        return hessian_.get();
      }

    protected:
      ScalarT           f_{};
      VectorT           x_;
      RealVectorT       x_lower_;
      RealVectorT       x_upper_;
      RealVectorT       gradient_;
      VectorT           g_;
      RealVectorT       g_lower_;
      RealVectorT       g_upper_;
      std::vector<bool> tag_;
      RealVectorT       lambda_;
      RealVectorT       z_lower_;
      RealVectorT       z_upper_;

      std::unique_ptr<CsrMatrixT> jacobian_;
      std::unique_ptr<CsrMatrixT> hessian_;
    };
  } // namespace Model
} // namespace GridKit
