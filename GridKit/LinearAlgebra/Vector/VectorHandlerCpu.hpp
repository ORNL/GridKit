#pragma once

#include <GridKit/MemoryUtilities/MemoryUtils.hpp>

namespace GridKit
{
  namespace LinearAlgebra
  {
    template <typename ScalarT, typename IdxT>
    class Vector;

    /**
     * @brief Implements dense (multi)vector kernels that operate on HOST data.
     *
     * All methods in this class assume the vectors passed in have up-to-date
     * data on HOST and operate directly on the raw HOST arrays.
     *
     * @author Slaven Peles <peless@ornl.gov>
     */
    template <typename ScalarT, typename IdxT>
    class VectorHandlerCpu
    {
    public:
      VectorHandlerCpu()  = default;
      ~VectorHandlerCpu() = default;

      // y = alpha * x + y
      void axpy(const ScalarT alpha, Vector<ScalarT, IdxT>* x, Vector<ScalarT, IdxT>* y);

      // dot: x \cdot y
      ScalarT dot(Vector<ScalarT, IdxT>* x, Vector<ScalarT, IdxT>* y);

      // scal: x = alpha * x
      void scal(const ScalarT alpha, Vector<ScalarT, IdxT>* x);

      // vector infinity norm
      ScalarT amax(Vector<ScalarT, IdxT>* x);

      // mass axpy: y = y - x*alpha, where x is [n x k] and alpha is [k x 1]; x is stored columnwise
      void axpyMulti(IdxT size, Vector<ScalarT, IdxT>* alpha, IdxT k, Vector<ScalarT, IdxT>* x, Vector<ScalarT, IdxT>* y);

      // mass dot: V^T x, where V is [n x k] and x is [k x 2], everything is stored and returned columnwise
      void dot2Multi(IdxT size, Vector<ScalarT, IdxT>* V, IdxT k, Vector<ScalarT, IdxT>* x, Vector<ScalarT, IdxT>* res);

      /** gemv:
       * if `transpose = N` (no), `x = beta*x +  alpha*V*y`,
       * where `x` is `[n x 1]`, `V` is `[n x k]` and `y` is `[k x 1]`.
       * if `transpose = T` (yes), `x = beta*x + alpha*V^T*y`,
       * where `x` is `[k x 1]`, `V` is `[n x k]` and `y` is `[n x 1]`.
       */
      void gemv(char                   transpose,
                IdxT                   k,
                const ScalarT          alpha,
                const ScalarT          beta,
                Vector<ScalarT, IdxT>* V,
                Vector<ScalarT, IdxT>* y,
                Vector<ScalarT, IdxT>* x);

      // Scale a vector by a diagonal matrix
      void scal(Vector<ScalarT, IdxT>* diag, Vector<ScalarT, IdxT>* vec);
      void scal(Vector<ScalarT, IdxT>* diag, Vector<ScalarT, IdxT>* vec, IdxT diag_offset);

      // Divide the elements of a vector by the elements of another vector
      int diagSolve(Vector<ScalarT, IdxT>* diag, Vector<ScalarT, IdxT>* vec);

      // Compute element-wise max of two vectors
      int max(Vector<ScalarT, IdxT>* x, Vector<ScalarT, IdxT>* y, Vector<ScalarT, IdxT>* out);

      // Compute element-wise absolute value of a vector
      int abs(Vector<ScalarT, IdxT>* in, Vector<ScalarT, IdxT>* out);
    };

  } // namespace LinearAlgebra
} // namespace GridKit
