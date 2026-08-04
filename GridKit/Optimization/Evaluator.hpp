#pragma once

#include <GridKit/LinearAlgebra/SparseMatrix/CsrMatrix.hpp>
#include <GridKit/LinearAlgebra/Vector/Vector.hpp>
#include <GridKit/ScalarTraits.hpp>

namespace GridKit
{
  namespace Optimization
  {
    /**
     * @brief Abstract nonlinear-program evaluator.
     *
     * Derivative structures are fixed by allocate() and remain unchanged until
     * the model is allocated again. For n variables and m constraints,
     * allocate() provides n-entry variable, bound, and objective-gradient
     * vectors; m-entry constraint and bound vectors; an m-by-n Jacobian; and an
     * n-by-n Hessian. Jacobian and Hessian values must correspond to the fixed
     * entries stored by their CSR matrices.
     */
    template <class scalar_type, typename index_type>
    class Evaluator
    {
    public:
      using ScalarT    = scalar_type;
      using IdxT       = index_type;
      using RealT      = typename ScalarTraits<ScalarT>::RealT;
      using VectorT    = LinearAlgebra::Vector<ScalarT, IdxT>;
      using CsrMatrixT = LinearAlgebra::CsrMatrix<RealT, IdxT>;

      virtual ~Evaluator() = default;

      virtual int allocate()   = 0;
      virtual int initialize() = 0;

      virtual IdxT size()            = 0;
      virtual IdxT sizeConstraints() = 0;

      virtual VectorT&       variables()       = 0;
      virtual const VectorT& variables() const = 0;

      virtual const VectorT& variableLowerBounds() const = 0;
      virtual const VectorT& variableUpperBounds() const = 0;

      virtual int   evaluateObjective() = 0;
      virtual RealT objective() const   = 0;

      virtual int            evaluateObjectiveGradient() = 0;
      virtual const VectorT& objectiveGradient() const   = 0;

      virtual int            evaluateConstraints()         = 0;
      virtual const VectorT& constraints() const           = 0;
      virtual const VectorT& constraintLowerBounds() const = 0;
      virtual const VectorT& constraintUpperBounds() const = 0;

      virtual int         evaluateJacobian()     = 0;
      virtual CsrMatrixT* getCsrJacobian() const = 0;

      /**
       * @brief Whether an exact, structurally invariant Jacobian is available.
       *
       * Returning false is a configuration error for the exact Ipopt adapter.
       * A model with only constant constraints returns true and may provide an
       * empty Jacobian structure.
       */
      virtual bool hasJacobian() = 0;

      /**
       * @brief Evaluate the exact lower-triangular Lagrangian Hessian.
       *
       * Computes
       * @f$ \sigma \nabla^2 f(x) + \sum_i \lambda_i \nabla^2 g_i(x) @f$.
       */
      virtual int evaluateHessian(RealT        objective_factor,
                                  const RealT* multipliers,
                                  IdxT         multiplier_count) = 0;

      virtual CsrMatrixT* getCsrHessian() const = 0;

      /**
       * @brief Whether an exact, structurally invariant Hessian is available.
       *
       * Returning false is a configuration error for the exact Ipopt adapter;
       * it never requests an approximate Hessian. A linear model returns true
       * and may provide an empty Hessian structure.
       */
      virtual bool hasHessian() = 0;
    };
  } // namespace Optimization
} // namespace GridKit
