#pragma once

#include <GridKit/LinearAlgebra/SparseMatrix/CsrMatrix.hpp>
#include <GridKit/ScalarTraits.hpp>

#include "GridKit/LinearAlgebra/Vector/Vector.hpp"

namespace GridKit
{
  namespace LinearAlgebra
  {

    /**
     * @brief Interface for linear solvers used by GridKit integrators, including
     *        \ref AnalysisManager::NativeDynamicSolver::Rosenbrock "Rosenbrock".
     *
     */
    template <class ScalarT, typename IdxT>
    class LinearSolver
    {
    public:
      using RealT = GridKit::ScalarTraits<ScalarT>::RealT;

      /**
       * @brief Configure the solver by analyzing the given matrix. The sparsity pattern of the given matrix must be set, and future calls
       * to other methods will use this sparsity pattern, as well as the data pointer of the given matrix.
       *
       * @param matrix The matrix to configure the solver with.
       * @return int An error code, or 0 if none.
       */
      virtual int configureSolver(GridKit::LinearAlgebra::CsrMatrix<RealT, IdxT>& matrix)                                       = 0;
      /**
       * @brief Setup the solver with matrix data.
       *
       * @param reuse_factors If true, factors from a previous call to `setupSolver()` are reused.
       * @pre \ref configureSolver() must be called first and every time the sparsity pattern of the matrix changes.
       * @pre If `reuse_factors` is true, then `setupSolver` must have been called once before with the same sparsity pattern.
       * @return int An error code, or 0 if none.
       */
      virtual int setupSolver(bool reuse_factors = false)                                                                       = 0;
      /**
       * @brief Perform a linear solve using the configured matrix.
       *
       * @param rhs The right-hand-side vector (input). No data is changed.
       * @param lhs The left-hand-side vector (output).
       * @pre \ref setupSolver must be called first and every time the data of the matrix changes.
       * @return int An error code, or 0 if none.
       *
       * @todo The data of the `rhs` parameter doesn't change, but we can't mark as cosntant because ReSolve currently requires
       * a non-`const` pointer for vectors, even if it will only read from that vector. In the future this might change, so this
       * can be updated to be `const`.
       */
      virtual int solve(GridKit::LinearAlgebra::Vector<ScalarT, IdxT>& rhs, GridKit::LinearAlgebra::Vector<ScalarT, IdxT>& lhs) = 0;
    };
  } // namespace LinearAlgebra
} // namespace GridKit
