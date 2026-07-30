#pragma once

#include <GridKit/LinearAlgebra/Vector/VectorHandlerCpu.hpp>
#include <GridKit/MemoryUtilities/MemoryUtils.hpp>

namespace GridKit
{
  namespace LinearAlgebra
  {
    template <typename scalar_type, typename index_type>
    class Vector;

    /**
     * @brief Dispatches dense (multi)vector operations to the implementation
     * appropriate for the requested memory space.
     *
     * Vector kernels operating on HOST memory are always available. Support
     * for DEVICE memory space is provided only when GridKit is built with GPU
     * support enabled and a device implementation for the requested kernel is
     * available.
     *
     * @author Slaven Peles <peless@ornl.gov>
     */
    template <typename scalar_type, typename index_type>
    class VectorHandler
    {
    public:
      using ScalarT = scalar_type;
      using IdxT    = index_type;

      VectorHandler()  = default;
      ~VectorHandler() = default;

      // y := alpha * x + y
      void axpy(const ScalarT          alpha,
                Vector<ScalarT, IdxT>* x,
                Vector<ScalarT, IdxT>* y,
                memory::MemorySpace    memspace);

      // Dot product of two vectors
      ScalarT dot(Vector<ScalarT, IdxT>* x, Vector<ScalarT, IdxT>* y, memory::MemorySpace memspace);

      // Scale vector by scalar
      void scal(const ScalarT alpha, Vector<ScalarT, IdxT>* x, memory::MemorySpace memspace);

      // Scale vector by diagonal matrix represented as a vector (i.e., vec = diag*vec)
      void scal(Vector<ScalarT, IdxT>* diag, Vector<ScalarT, IdxT>* vec, memory::MemorySpace memspace);
      void scal(Vector<ScalarT, IdxT>* diag,
                Vector<ScalarT, IdxT>* vec,
                IdxT                   diag_offset,
                memory::MemorySpace    memspace);

      // axpy for multivectors
      void axpyMulti(IdxT                   size,
                     Vector<ScalarT, IdxT>* alpha,
                     IdxT                   k,
                     Vector<ScalarT, IdxT>* x,
                     Vector<ScalarT, IdxT>* y,
                     memory::MemorySpace    memspace);

      // Dot product of two vectors with a multivector V
      void dot2Multi(IdxT                   size,
                     Vector<ScalarT, IdxT>* V,
                     IdxT                   k,
                     Vector<ScalarT, IdxT>* x,
                     Vector<ScalarT, IdxT>* res,
                     memory::MemorySpace    memspace);

      // Dense matrix-vector product.
      void gemv(char                   transpose,
                IdxT                   k, // number of vectors from multivector V to use
                const ScalarT          alpha,
                const ScalarT          beta,
                Vector<ScalarT, IdxT>* V,
                Vector<ScalarT, IdxT>* y,
                Vector<ScalarT, IdxT>* x,
                memory::MemorySpace    memspace);

      int diagSolve(Vector<ScalarT, IdxT>* diag, Vector<ScalarT, IdxT>* vec, memory::MemorySpace memspace);
      int max(Vector<ScalarT, IdxT>* x, Vector<ScalarT, IdxT>* y, Vector<ScalarT, IdxT>* out, memory::MemorySpace memspace);
      int abs(Vector<ScalarT, IdxT>* in, Vector<ScalarT, IdxT>* out, memory::MemorySpace memspace);

      // Vector infinity norm
      ScalarT amax(Vector<ScalarT, IdxT>* x, memory::MemorySpace memspace);

    private:
      VectorHandlerCpu<ScalarT, IdxT> cpuImpl_;
    }; // class VectorHandler

  } // namespace LinearAlgebra
} // namespace GridKit
