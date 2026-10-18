/**
 * @file OptimizationEvaluator.hpp
 * @brief Abstract nonlinear program with exact sparse derivatives.
 */

#pragma once

#include <memory>

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
     * Minimize \f$f(x)\f$ subject to \f$c(x) = 0\f$,
     * \f$d_l \le d(x) \le d_u\f$, and \f$x_l \le x \le x_u\f$.
     * The Lagrangian is
     * \f$L = \sigma f + \lambda_c^T c + \lambda_d^T d\f$.
     * Jacobians and the lower Lagrangian Hessian retain the sparsity
     * patterns set by `allocate()`. Evaluations are at `x()`, and methods
     * returning `int` return 0 on success.
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

      /// Number the variables and constraints, and set the bounds, the starting point, and the sparsity patterns
      virtual int allocate() = 0;

      /// Set the starting point `x()`
      virtual int initialize() = 0;

      /// Evaluate \f$f(x)\f$ into `f()`
      virtual int evaluateObjective() = 0;

      /// Evaluate \f$\nabla f(x)\f$ into `gradient()`
      virtual int evaluateGradient() = 0;

      /// Evaluate \f$c(x)\f$ into `c()`
      virtual int evaluateEqualities() = 0;

      /// Evaluate \f$d(x)\f$ into `d()`
      virtual int evaluateInequalities() = 0;

      /// Evaluate the Jacobian of \f$c(x)\f$ into `getCsrEqualityJacobian()`
      virtual int evaluateEqualityJacobian() = 0;

      /// Evaluate the Jacobian of \f$d(x)\f$ into `getCsrInequalityJacobian()`
      virtual int evaluateInequalityJacobian() = 0;

      /**
       * @brief Evaluate the lower triangle of the Hessian of \f$L\f$ into
       * `getCsrHessian()`
       *
       * @param[in] sigma - Objective factor \f$\sigma\f$
       * @param[in] lambda_c - Equality multipliers \f$\lambda_c\f$, `sizeEqualities()` values
       * @param[in] lambda_d - Inequality multipliers \f$\lambda_d\f$, `sizeInequalities()` values
       */
      virtual int evaluateHessian(RealT        sigma,
                                  const RealT* lambda_c,
                                  const RealT* lambda_d) = 0;

      /// Number of variables
      IdxT size() const
      {
        return x_.getSize();
      }

      /// Number of equality constraints
      IdxT sizeEqualities() const
      {
        return c_.getSize();
      }

      /// Number of inequality constraints
      IdxT sizeInequalities() const
      {
        return d_.getSize();
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

      /// Equality constraints at the last `evaluateEqualities()`
      const VectorT& c() const
      {
        return c_;
      }

      /// Inequality constraints at the last `evaluateInequalities()`
      const VectorT& d() const
      {
        return d_;
      }

      /// Inequality constraint lower bounds
      const RealVectorT& dLower() const
      {
        return d_lower_;
      }

      /// Inequality constraint upper bounds
      const RealVectorT& dUpper() const
      {
        return d_upper_;
      }

      /// Jacobian of the equality constraints
      CsrMatrixT* getCsrEqualityJacobian() const
      {
        return equality_jacobian_.get();
      }

      /// Jacobian of the inequality constraints
      CsrMatrixT* getCsrInequalityJacobian() const
      {
        return inequality_jacobian_.get();
      }

      /// Lower triangle of the Lagrangian Hessian
      CsrMatrixT* getCsrHessian() const
      {
        return hessian_.get();
      }

    protected:
      ScalarT     f_{};
      VectorT     x_;
      RealVectorT x_lower_;
      RealVectorT x_upper_;
      RealVectorT gradient_;
      VectorT     c_;
      VectorT     d_;
      RealVectorT d_lower_;
      RealVectorT d_upper_;

      std::unique_ptr<CsrMatrixT> equality_jacobian_;
      std::unique_ptr<CsrMatrixT> inequality_jacobian_;
      std::unique_ptr<CsrMatrixT> hessian_;
    };
  } // namespace Model
} // namespace GridKit
