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
     * Minimize \f$f(x)\f$ subject to \f$g_l \le g(x) \le g_u\f$ and
     * \f$x_l \le x \le x_u\f$. The constraint Jacobian and the lower triangle
     * of the Lagrangian Hessian keep the sparsity pattern set by `allocate()`.
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

      virtual int allocate()   = 0;
      virtual int initialize() = 0;

      virtual int evaluateObjective()   = 0;
      virtual int evaluateGradient()    = 0;
      virtual int evaluateConstraints() = 0;
      virtual int evaluateJacobian()    = 0;

      /// Evaluate the Hessian of \f$\sigma f(x) + \lambda^T g(x)\f$
      virtual int evaluateHessian(RealT sigma, const RealT* lambda) = 0;

      IdxT size() const
      {
        return x_.getSize();
      }

      IdxT sizeConstraints() const
      {
        return g_.getSize();
      }

      const ScalarT& objective() const
      {
        return f_;
      }

      VectorT& x()
      {
        return x_;
      }

      const VectorT& x() const
      {
        return x_;
      }

      const RealVectorT& xLower() const
      {
        return x_lower_;
      }

      const RealVectorT& xUpper() const
      {
        return x_upper_;
      }

      const RealVectorT& gradient() const
      {
        return gradient_;
      }

      const VectorT& g() const
      {
        return g_;
      }

      const RealVectorT& gLower() const
      {
        return g_lower_;
      }

      const RealVectorT& gUpper() const
      {
        return g_upper_;
      }

      CsrMatrixT* getCsrJacobian() const
      {
        return jacobian_.get();
      }

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
      VectorT     g_;
      RealVectorT g_lower_;
      RealVectorT g_upper_;

      std::unique_ptr<CsrMatrixT> jacobian_;
      std::unique_ptr<CsrMatrixT> hessian_;
    };
  } // namespace Model
} // namespace GridKit
