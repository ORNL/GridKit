#pragma once

#include <memory>

#include <GridKit/LinearAlgebra/Solver/LinearSolver.hpp>

#include <resolve/MemoryUtils.hpp>
#include <resolve/SystemSolver.hpp>
#include <resolve/matrix/Csr.hpp>

namespace GridKit
{
  namespace LinearAlgebra
  {
    template <typename scalar_type, typename index_type>
    class ResolveSystemSolver : public LinearSolver<scalar_type, index_type>
    {
    public:
      using ScalarT = scalar_type;
      using IdxT    = index_type;
      using RealT   = LinearSolver<ScalarT, IdxT>::RealT;

      ResolveSystemSolver(ReSolve::SystemSolver& lin_solver, GridKit::memory::MemorySpace memspace = GridKit::memory::HOST);

      int configureSolver(GridKit::LinearAlgebra::CsrMatrix<RealT, IdxT>& matrix) final;
      int setupSolver(bool reuse_factors = false) final;
      int solve(GridKit::LinearAlgebra::Vector<ScalarT, IdxT>& rhs, GridKit::LinearAlgebra::Vector<ScalarT, IdxT>& lhs) final;

    private:
      /**
       * @brief The ReSolve matrix used for the linear solve. Its data pointers should just point to a \ref GridKit::LinearAlgebra::CsrMatrix.
       *
       */
      std::unique_ptr<ReSolve::matrix::Csr> matrix_;
      /**
       * @brief The ReSolve linear solver to use.
       *
       */
      ReSolve::SystemSolver&                lin_solver_;
      /**
       * @brief The memspace the data is stored in.
       *
       */
      ReSolve::memory::MemorySpace          memspace_;
    };
  } // namespace LinearAlgebra
} // namespace GridKit
