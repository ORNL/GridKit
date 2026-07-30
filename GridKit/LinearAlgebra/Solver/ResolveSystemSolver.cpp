#include "ResolveSystemSolver.hpp"

#include "GridKit/MemoryUtilities/MemoryUtils.hpp"
#include <resolve/vector/Vector.hpp>

namespace GridKit
{
  namespace LinearAlgebra
  {
    /**
     * @brief Helper function to convert GridKit `MemorySpace` into ReSolve `MemorySpace`.
     *
     */
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

    template <typename scalar_type, typename index_type>
    ResolveSystemSolver<scalar_type, index_type>::ResolveSystemSolver(ReSolve::SystemSolver& lin_solver, GridKit::memory::MemorySpace memspace)
      : lin_solver_(lin_solver), memspace_(memorySpaceAsResolve(memspace))
    {
    }

    /**
     * @brief Creates a new ReSolve matrix with the dimensions of the given matrix, then sets its data pointers to those of the given matrix.
     * Sets up the linear solver to use this new ReSolve matrix.
     *
     * @todo Right now preconditioning doesn't work. There should be a ReSolve PR soon for this.
     *
     */
    template <typename scalar_type, typename index_type>
    int ResolveSystemSolver<scalar_type, index_type>::configureSolver(GridKit::LinearAlgebra::CsrMatrix<RealT, IdxT>& matrix)
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

    template <typename scalar_type, typename index_type>
    int ResolveSystemSolver<scalar_type, index_type>::setupSolver(bool reuse_factors)
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

    /**
     * @brief Perform the solver by creating two new ReSolve vectors and setting their data to point to their GridKit counterparts.
     * Then ask ReSolve to perform the solve. Since the matrix should be pointing at the correct GridKit buffers and the vectors
     * as well, this solve should correctly fill `lhs`'s data buffer.
     *
     */
    template <typename scalar_type, typename index_type>
    int ResolveSystemSolver<scalar_type, index_type>::solve(GridKit::LinearAlgebra::Vector<ScalarT, IdxT>& rhs, GridKit::LinearAlgebra::Vector<ScalarT, IdxT>& lhs)
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
