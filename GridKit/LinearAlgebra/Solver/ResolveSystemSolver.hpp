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
    template <class ScalarT, typename IdxT>
    class ResolveSystemSolver : public LinearSolver<ScalarT, IdxT>
    {
    public:
      using RealT = LinearSolver<ScalarT, IdxT>::RealT;

      ResolveSystemSolver(ReSolve::SystemSolver& lin_solver, GridKit::memory::MemorySpace memspace = GridKit::memory::HOST);

      int configureSolver(GridKit::LinearAlgebra::CsrMatrix<RealT, IdxT>& matrix) final;
      int setupSolver(bool reuse_factors = false) final;
      int solve(GridKit::LinearAlgebra::Vector<ScalarT, IdxT>& rhs, GridKit::LinearAlgebra::Vector<ScalarT, IdxT>& lhs) final;

    private:
      std::unique_ptr<ReSolve::matrix::Csr> matrix_;
      ReSolve::SystemSolver&                lin_solver_;
      ReSolve::memory::MemorySpace          memspace_;
    };
  } // namespace LinearAlgebra
} // namespace GridKit