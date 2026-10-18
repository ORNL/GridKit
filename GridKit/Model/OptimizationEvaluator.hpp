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
     * @brief Abstract class describing a nonlinear program
     *
     * Minimize \f$f(x)\f$ subject to \f$g_l \le g(x) \le g_u\f$ and
     * \f$x_l \le x \le x_u\f$. `tag()` marks the inequality constraints; the
     * others are equalities with \f$g_l = g_u = 0\f$. The Lagrangian is
     * \f$L = \sigma f + \lambda^T g\f$. `allocate()` sets the sizes, the tags,
     * and the sparsity patterns, which the Jacobian and the lower Lagrangian
     * Hessian retain. `initialize()` sets the bounds and the initial values.
     * Evaluations are at `x()` and `lambda()`, and methods returning `int`
     * return 0 on success.
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

      /// Number the variables and constraints, tag the inequalities, and set the sparsity patterns
      virtual int allocate() = 0;

      /// Fill `tag()`, true for each inequality constraint
      virtual int tagInequality() = 0;

      /// Set the bounds, the initial values of the variables `x()`, and zero multipliers
      virtual int initialize() = 0;

      /// Evaluate \f$f(x)\f$ into `f()`
      virtual int evaluateObjective() = 0;

      /// Evaluate \f$\nabla f(x)\f$ into `gradient()`
      virtual int evaluateGradient() = 0;

      /// Evaluate \f$g(x)\f$ into `g()`
      virtual int evaluateConstraints() = 0;

      /// Evaluate the Jacobian of \f$g(x)\f$ into `getCsrJacobian()`
      virtual int evaluateJacobian() = 0;

      /**
       * @brief Evaluate the lower triangle of the Hessian of \f$L\f$ at the
       * multipliers `lambda()` into `getCsrHessian()`
       *
       * @param[in] sigma - Objective factor \f$\sigma\f$
       */
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

      /// Objective at the last `evaluateObjective()`
      const ScalarT& f() const
      {
        return f_;
      }

      /// Variables
      VectorT& x()
      {
        return x_;
      }

      /// Variables
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

      /// Objective gradient at the last `evaluateGradient()`
      const RealVectorT& gradient() const
      {
        return gradient_;
      }

      /// Constraints at the last `evaluateConstraints()`
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

      /// Inequality tag of each constraint, set by `tagInequality()`
      const std::vector<bool>& tag() const
      {
        return tag_;
      }

      /// Constraint multipliers \f$\lambda\f$
      RealVectorT& lambda()
      {
        return lambda_;
      }

      /// Constraint multipliers \f$\lambda\f$
      const RealVectorT& lambda() const
      {
        return lambda_;
      }

      /// Multipliers of the variable lower bounds
      RealVectorT& zLower()
      {
        return z_lower_;
      }

      /// Multipliers of the variable lower bounds
      const RealVectorT& zLower() const
      {
        return z_lower_;
      }

      /// Multipliers of the variable upper bounds
      RealVectorT& zUpper()
      {
        return z_upper_;
      }

      /// Multipliers of the variable upper bounds
      const RealVectorT& zUpper() const
      {
        return z_upper_;
      }

      /// Jacobian of the constraints
      CsrMatrixT* getCsrJacobian() const
      {
        return jacobian_.get();
      }

      /// Lower triangle of the Lagrangian Hessian
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
