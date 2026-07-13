#include "ResolveSystemSolver.hpp"

#include "GridKit/MemoryUtilities/MemoryUtils.hpp"
#include <resolve/vector/Vector.hpp>

namespace GridKit
{
  namespace LinearAlgebra
  {
    ReSolve::memory::MemorySpace memorySpaceAsResolve(memory::MemorySpace memspace)
    {
      switch (memspace)
      {
      case memory::HOST:
        return ReSolve::memory::HOST;
      case memory::DEVICE:
        return ReSolve::memory::DEVICE;
      default:
        throw "Memory space not supported";
      }
    }

    template <class ScalarT, typename IdxT>
    ResolveSystemSolver<ScalarT, IdxT>::ResolveSystemSolver(ReSolve::SystemSolver& lin_solver, GridKit::memory::MemorySpace memspace) : lin_solver_(lin_solver), memspace_(memorySpaceAsResolve(memspace))
    {
    }

    template <class ScalarT, typename IdxT>
    int ResolveSystemSolver<ScalarT, IdxT>::configureSolver(GridKit::LinearAlgebra::CsrMatrix<RealT, IdxT>& matrix)
    {
      matrix_ = std::make_unique<ReSolve::matrix::Csr>(matrix.getNumRows(), matrix.getNumColumns(), matrix.getNnz());

      if (int err = matrix_->setDataPointers(matrix.getRowData(), matrix.getColData(), matrix.getValues(), memspace_))
        return err;

      if (int err = lin_solver_.setMatrix(matrix_.get()))
        return err;

      if (int err = lin_solver_.analyze())
        return err;

      // TODO: This function needs to be called to properly use a preconditioner in ReSolve (if there is any), but currently will error
      // if there is no preconditioner configured. Once we can detect if a preconditioner is configured, we can restore this functionality.
      // Also, we will always want to use *right* preconditioning.
      // BUBBLE_FAIL(lin_solver_.preconditionerSetup("right"));

      return 0;
    }

    template <class ScalarT, typename IdxT>
    int ResolveSystemSolver<ScalarT, IdxT>::setupSolver(bool reuse_factors)
    {
      if (reuse_factors)
      {
        if (int err = lin_solver_.refactorize())
          return err;
      }
      else
      {
        if (int err = lin_solver_.factorize())
          return err;
      }

      return 0;
    }

    template <class ScalarT, typename IdxT>
    int ResolveSystemSolver<ScalarT, IdxT>::solve(GridKit::LinearAlgebra::Vector<ScalarT, IdxT>& rhs, GridKit::LinearAlgebra::Vector<ScalarT, IdxT>& lhs)
    {
      ReSolve::vector::Vector resolve_rhs(rhs.getSize());
      ReSolve::vector::Vector resolve_lhs(lhs.getSize());

      if (int err = resolve_rhs.setData(rhs.getData(), memspace_))
        return err;

      if (int err = resolve_lhs.setData(lhs.getData(), memspace_))
        return err;

      if (int err = lin_solver_.solve(&resolve_rhs, &resolve_lhs))
        return err;

      return 0;
    }

    template class ResolveSystemSolver<double, int>;
  } // namespace LinearAlgebra
} // namespace GridKit