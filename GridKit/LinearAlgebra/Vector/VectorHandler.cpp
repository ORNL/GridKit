#include "VectorHandler.hpp"

#include <cassert>
#include <cmath>

#include <GridKit/LinearAlgebra/Vector/Vector.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace LinearAlgebra
  {
    using out = GridKit::Utilities::Logger;

    /**
     * @brief dot product of two vectors i.e, a = x^Ty
     *
     * @param[in] x The first vector
     * @param[in] y The second vector
     * @param[in] memspace Memory space the operation is computed in (HOST or DEVICE)
     *
     * @return dot product of _x_ and _y_
     */
    template <typename ScalarT, typename IdxT>
    ScalarT VectorHandler<ScalarT, IdxT>::dot(Vector<ScalarT, IdxT>* x, Vector<ScalarT, IdxT>* y, memory::MemorySpace memspace)
    {
      switch (memspace)
      {
      case memory::HOST:
        return cpuImpl_.dot(x, y);
      case memory::DEVICE:
        out::error() << "VectorHandler::dot - DEVICE memory space not yet supported\n";
        return static_cast<ScalarT>(NAN);
      }
      return static_cast<ScalarT>(NAN);
    }

    /**
     * @brief scale a vector by a constant i.e, x = alpha*x where alpha is a constant
     *
     * @param[in] alpha The constant
     * @param[in,out] x The vector
     * @param[in] memspace Memory space the operation is computed in (HOST or DEVICE)
     */
    template <typename ScalarT, typename IdxT>
    void VectorHandler<ScalarT, IdxT>::scal(const ScalarT alpha, Vector<ScalarT, IdxT>* x, memory::MemorySpace memspace)
    {
      switch (memspace)
      {
      case memory::HOST:
        cpuImpl_.scal(alpha, x);
        break;
      case memory::DEVICE:
        out::error() << "VectorHandler::scal - DEVICE memory space not yet supported\n";
        break;
      }
    }

    /**
     * @brief compute infinity norm of a vector (i.e., find an entry with largest absolute value)
     *
     * @param[in] x The vector
     * @param[in] memspace Memory space the operation is computed in (HOST or DEVICE)
     *
     * @return infinity norm of _x_
     */
    template <typename ScalarT, typename IdxT>
    ScalarT VectorHandler<ScalarT, IdxT>::amax(Vector<ScalarT, IdxT>* x, memory::MemorySpace memspace)
    {
      switch (memspace)
      {
      case memory::HOST:
        return cpuImpl_.amax(x);
      case memory::DEVICE:
        out::error() << "VectorHandler::amax - DEVICE memory space not yet supported\n";
        return static_cast<ScalarT>(NAN);
      }
      return static_cast<ScalarT>(NAN);
    }

    /**
     * @brief axpy i.e, y = alpha*x+y where alpha is a constant
     *
     * @param[in] alpha The constant
     * @param[in] x The first vector
     * @param[in,out] y The second vector (result is returned in y)
     * @param[in] memspace Memory space the operation is computed in (HOST or DEVICE)
     */
    template <typename ScalarT, typename IdxT>
    void VectorHandler<ScalarT, IdxT>::axpy(const ScalarT          alpha,
                                            Vector<ScalarT, IdxT>* x,
                                            Vector<ScalarT, IdxT>* y,
                                            memory::MemorySpace    memspace)
    {
      switch (memspace)
      {
      case memory::HOST:
        cpuImpl_.axpy(alpha, x, y);
        break;
      case memory::DEVICE:
        out::error() << "VectorHandler::axpy - DEVICE memory space not yet supported\n";
        break;
      }
    }

    /**
     * @brief gemv computes dense matrix-vector product.
     *
     * If `transpose = N` (no), `x := beta*x + alpha*V*y`,
     * where `x` is `[n x 1]`, `V` is `[n x k]` and `y` is `[k x 1]`.
     * If `transpose = T` (yes), `x := beta*x + alpha*V^T*y`,
     * where `x` is `[k x 1]`, `V` is `[n x k]` and `y` is `[n x 1]`.
     *
     * @param[in] transpose - yes (T) or no (N)
     * @param[in] k         - Number of columns in (non-transposed) matrix to use
     * @param[in] alpha     - Constant scalar
     * @param[in] beta      - Constant scalar
     * @param[in] V         - Multivector containing the matrix, organized columnwise
     * @param[in] y         - Vector, k x 1 if N and n x 1 if T
     * @param[in,out] x     - Vector, n x 1 if N and k x 1 if T
     * @param[in] memspace  - Memory space the operation is computed in (HOST or DEVICE)
     *
     * @note Parameter k is not the total number of columns in V but the number
     * of columns to use in matrix-vector product.
     */
    template <typename ScalarT, typename IdxT>
    void VectorHandler<ScalarT, IdxT>::gemv(char                   transpose,
                                            IdxT                   k,
                                            const ScalarT          alpha,
                                            const ScalarT          beta,
                                            Vector<ScalarT, IdxT>* V,
                                            Vector<ScalarT, IdxT>* y,
                                            Vector<ScalarT, IdxT>* x,
                                            memory::MemorySpace    memspace)
    {
      switch (memspace)
      {
      case memory::HOST:
        cpuImpl_.gemv(transpose, k, alpha, beta, V, y, x);
        break;
      case memory::DEVICE:
        out::error() << "VectorHandler::gemv - DEVICE memory space not yet supported\n";
        break;
      }
    }

    /**
     * @brief Multivector axpy: y := y - \sum_i alpha_i x_i
     *
     * @param[in] size number of elements in y
     * @param[in] alpha vector size k x 1
     * @param[in] k number of vectors in the multivector x
     * @param[in] x (multi)vector [size x k]
     * @param[in,out] y vector size size x 1 (this is where the result is stored)
     * @param[in] memspace Memory space the operation is computed in (HOST or DEVICE)
     *
     * @pre _k_ > 0, _size_ > 0, _size_ = x->getSize()
     */
    template <typename ScalarT, typename IdxT>
    void VectorHandler<ScalarT, IdxT>::axpyMulti(IdxT                   size,
                                                 Vector<ScalarT, IdxT>* alpha,
                                                 IdxT                   k,
                                                 Vector<ScalarT, IdxT>* x,
                                                 Vector<ScalarT, IdxT>* y,
                                                 memory::MemorySpace    memspace)
    {
      assert(y->getSize() == x->getSize() && "Sizes of x and y must match!\n");
      assert(alpha->getSize() == k && "Size of alpha must match k!\n");

      switch (memspace)
      {
      case memory::HOST:
        cpuImpl_.axpyMulti(size, alpha, k, x, y);
        break;
      case memory::DEVICE:
        out::error() << "VectorHandler::axpyMulti - DEVICE memory space not yet supported\n";
        break;
      }
    }

    /**
     * @brief Multivector dot product, i.e V^T x
     *
     * Computes V^T x with k vectors from multivector V. Result is stored in `res`.
     *
     * @param[in] size     - Number of elements in a single vector in V
     * @param[in] V        - Multivector; k vectors of size n x 1 each
     * @param[in] k        - Number of vectors in V to use
     * @param[in] x        - Multivector; 2 vectors of size n x 1 each
     * @param[out] res     - Multivector; 2 vectors size k x 1 each
     * @param[in] memspace - Memory space the operation is computed in (HOST or DEVICE)
     *
     * @pre _size_ > 0, _k_ > 0, size = x->getSize().
     * @pre _res_ needs to be allocated to k x 2 size.
     */
    template <typename ScalarT, typename IdxT>
    void VectorHandler<ScalarT, IdxT>::dot2Multi(IdxT                   size,
                                                 Vector<ScalarT, IdxT>* V,
                                                 IdxT                   k,
                                                 Vector<ScalarT, IdxT>* x,
                                                 Vector<ScalarT, IdxT>* res,
                                                 memory::MemorySpace    memspace)
    {
      assert(x->getSize() == V->getSize() && "Sizes of V and x do not match!\n");
      assert(res->getSize() == k && "Size of `res` must match k!\n");

      switch (memspace)
      {
      case memory::HOST:
        cpuImpl_.dot2Multi(size, V, k, x, res);
        break;
      case memory::DEVICE:
        out::error() << "VectorHandler::dot2Multi - DEVICE memory space not yet supported\n";
        break;
      }
    }

    /**
     * @brief Scale a vector by a diagonal matrix
     *
     * @param[in] diag - vector representing the diagonal matrix
     * @param[in,out] vec - vector to be scaled
     * @param[in] memspace - Memory space the operation is computed in (HOST or DEVICE)
     *
     * @pre The diagonal vector must be of the same size as the vector.
     */
    template <typename ScalarT, typename IdxT>
    void VectorHandler<ScalarT, IdxT>::scal(Vector<ScalarT, IdxT>* diag, Vector<ScalarT, IdxT>* vec, memory::MemorySpace memspace)
    {
      assert(diag->getSize() == vec->getSize() && "Diagonal vector must be of the same size as the vector.");

      switch (memspace)
      {
      case memory::HOST:
        cpuImpl_.scal(diag, vec);
        break;
      case memory::DEVICE:
        out::error() << "VectorHandler::scal - DEVICE memory space not yet supported\n";
        break;
      }
    }

    /**
     * @brief Scale a vector by a diagonal matrix represented by a contiguous subvector of an input vector
     *
     * @param[in] diag - vector representing the diagonal matrix
     * @param[in,out] vec - vector to be scaled
     * @param[in] diag_offset - the index of diag where the diagonal matrix begins (inclusive)
     * @param[in] memspace - Memory space the operation is computed in (HOST or DEVICE)
     */
    template <typename ScalarT, typename IdxT>
    void VectorHandler<ScalarT, IdxT>::scal(Vector<ScalarT, IdxT>* diag,
                                            Vector<ScalarT, IdxT>* vec,
                                            IdxT                   diag_offset,
                                            memory::MemorySpace    memspace)
    {
      switch (memspace)
      {
      case memory::HOST:
        cpuImpl_.scal(diag, vec, diag_offset);
        break;
      case memory::DEVICE:
        out::error() << "VectorHandler::scal - DEVICE memory space not yet supported\n";
        break;
      }
    }

    /**
     * @brief Multiplies vector by an inverse of a diagonal matrix.
     *
     * @param[in]  diag   - diagonal matrix stored in a vector object
     * @param[in,out] vec - vector to be divided
     * @param[in] memspace - Memory space the operation is computed in (HOST or DEVICE)
     *
     * @pre The two vectors must be the same size
     *
     * @return 0 if successful, 1 otherwise
     */
    template <typename ScalarT, typename IdxT>
    int VectorHandler<ScalarT, IdxT>::diagSolve(Vector<ScalarT, IdxT>* diag, Vector<ScalarT, IdxT>* vec, memory::MemorySpace memspace)
    {
      assert(diag->getSize() == vec->getSize() && "Diagonal vector must be of the same size as the vector.");

      switch (memspace)
      {
      case memory::HOST:
        return cpuImpl_.diagSolve(diag, vec);
      case memory::DEVICE:
        out::error() << "VectorHandler::diagSolve - DEVICE memory space not yet supported\n";
        return 1;
      }
      return 1;
    }

    /**
     * @brief Takes the element-wise max between two vectors.
     *
     * @param[in]  x        - The first vector
     * @param[in]  y        - The second vector
     * @param[out] out      - The output vector
     * @param[in]  memspace - Memory space the operation is computed in (HOST or DEVICE)
     *
     * @pre The two vectors must be the same size
     *
     * @return 0 if successful, 1 otherwise
     */
    template <typename ScalarT, typename IdxT>
    int VectorHandler<ScalarT, IdxT>::max(Vector<ScalarT, IdxT>* x, Vector<ScalarT, IdxT>* y, Vector<ScalarT, IdxT>* out, memory::MemorySpace memspace)
    {
      assert(x->getSize() == y->getSize() && "Vectors must be the same size.");
      assert(x->getSize() == out->getSize() && "Vectors must be the same size.");

      switch (memspace)
      {
      case memory::HOST:
        return cpuImpl_.max(x, y, out);
      case memory::DEVICE:
        GridKit::LinearAlgebra::out::error() << "VectorHandler::max - DEVICE memory space not yet supported\n";
        return 1;
      }
      return 1;
    }

    /**
     * @brief Computes the element-wise absolute value of a vector.
     *
     * @param[in]  in       - Input vector
     * @param[out] out      - Output vector
     * @param[in]  memspace - Memory space the operation is computed in (HOST or DEVICE)
     *
     * @return 0 if successful, 1 otherwise
     */
    template <typename ScalarT, typename IdxT>
    int VectorHandler<ScalarT, IdxT>::abs(Vector<ScalarT, IdxT>* in, Vector<ScalarT, IdxT>* out, memory::MemorySpace memspace)
    {
      assert(in->getSize() == out->getSize() && "Vector sizes do not match!\n");

      switch (memspace)
      {
      case memory::HOST:
        return cpuImpl_.abs(in, out);
      case memory::DEVICE:
        GridKit::LinearAlgebra::out::error() << "VectorHandler::abs - DEVICE memory space not yet supported\n";
        return 1;
      }
      return 1;
    }

    // Available template instantiations
    template class VectorHandler<double, size_t>;
    template class VectorHandler<double, int>;

  } // namespace LinearAlgebra
} // namespace GridKit
