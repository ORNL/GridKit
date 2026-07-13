#pragma once

#include <GridKit/LinearAlgebra/SparseMatrix/CsrMatrix.hpp>
#include <GridKit/ScalarTraits.hpp>

#include "GridKit/LinearAlgebra/Vector/Vector.hpp"

namespace GridKit
{
  namespace LinearAlgebra
  {

    template <class ScalarT, typename IdxT>
    class LinearSolver
    {
    public:
      using RealT = GridKit::ScalarTraits<ScalarT>::RealT;

      virtual int configureSolver(GridKit::LinearAlgebra::CsrMatrix<RealT, IdxT>& matrix)                                       = 0;
      virtual int setupSolver(bool reuse_factors = false)                                                                       = 0;
      virtual int solve(GridKit::LinearAlgebra::Vector<ScalarT, IdxT>& rhs, GridKit::LinearAlgebra::Vector<ScalarT, IdxT>& lhs) = 0;
    };
  } // namespace LinearAlgebra
} // namespace GridKit
