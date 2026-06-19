#include "VectorHandlerCpu.hpp"

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
     *
     * @return dot product of _x_ and _y_
     */
    template <typename scalar_type, typename index_type>
    scalar_type VectorHandlerCpu<scalar_type, index_type>::dot(Vector<ScalarT, IdxT>* x, Vector<ScalarT, IdxT>* y)
    {
      const ScalarT* x_data = x->getData(memory::HOST);
      const ScalarT* y_data = y->getData(memory::HOST);
      ScalarT        sum    = 0.0;
      ScalarT        c      = 0.0;
      // Kahan summation to reduce round-off error accumulation
      for (IdxT i = 0; i < x->getSize(); ++i)
      {
        ScalarT y = (x_data[i] * y_data[i]) - c;
        ScalarT t = sum + y;
        c         = (t - sum) - y;
        sum       = t;
      }
      return sum;
    }

    /**
     * @brief scale a vector by a constant i.e, x = alpha*x where alpha is a constant
     *
     * @param[in] alpha The constant
     * @param[in,out] x The vector
     */
    template <typename scalar_type, typename index_type>
    void VectorHandlerCpu<scalar_type, index_type>::scal(const ScalarT alpha, Vector<ScalarT, IdxT>* x)
    {
      ScalarT* x_data = x->getData(memory::HOST);

      for (IdxT i = 0; i < x->getSize(); ++i)
      {
        x_data[i] *= alpha;
      }
      x->setDataUpdated(memory::HOST);
    }

    /**
     * @brief compute infinity norm of a vector (i.e., find an entry with largest absolute value)
     *
     * @param[in] x vector
     *
     * @return infinity norm of _x_
     */
    template <typename scalar_type, typename index_type>
    scalar_type VectorHandlerCpu<scalar_type, index_type>::amax(Vector<ScalarT, IdxT>* x)
    {
      const ScalarT* x_data = x->getData(memory::HOST);

      ScalarT vecmax = std::abs(x_data[0]);
      ScalarT v;
      for (IdxT i = 1; i < x->getSize(); ++i)
      {
        v = std::abs(x_data[i]);
        if (v > vecmax)
        {
          vecmax = v;
        }
      }
      return vecmax;
    }

    /**
     * @brief axpy i.e, y = alpha*x+y where alpha is a constant
     *
     * @param[in] alpha The constant
     * @param[in] x The first vector
     * @param[in,out] y The second vector (result is returned in y)
     */
    template <typename scalar_type, typename index_type>
    void VectorHandlerCpu<scalar_type, index_type>::axpy(const ScalarT alpha, Vector<ScalarT, IdxT>* x, Vector<ScalarT, IdxT>* y)
    {
      ScalarT* x_data = x->getData(memory::HOST);
      ScalarT* y_data = y->getData(memory::HOST);
      for (IdxT i = 0; i < x->getSize(); ++i)
      {
        y_data[i] = alpha * x_data[i] + y_data[i];
      }
      y->setDataUpdated(memory::HOST);
    }

    /**
     * @brief gemv computes matrix-vector product where the matrix and vectors are dense.
     *        i.e., x = beta*x + alpha*V*y
     *
     * @param[in] transpose - transposed = 'T' or not 'N'
     * @param[in] k Number of columns in (non-transposed) matrix to use
     * @param[in] alpha Constant scalar
     * @param[in] beta Constant scalar
     * @param[in] V Multivector containing the matrix, organized columnwise
     * @param[in] y Vector, k x 1 if N and n x 1 if T
     * @param[in,out] x Vector, n x 1 if N and k x 1 if T
     *
     * @note Parameter k is not the total number of columns in V but the number
     * of columns to use in matrix-vector product.
     *
     * @pre Number of columns in V >= k
     * @pre If transpose = N, size of y must equal k. If transpose = T, size of
     * x must equal k.
     */
    template <typename scalar_type, typename index_type>
    void VectorHandlerCpu<scalar_type, index_type>::gemv(char                   transpose,
                                                         IdxT                   k,
                                                         const ScalarT          alpha,
                                                         const ScalarT          beta,
                                                         Vector<ScalarT, IdxT>* V,
                                                         Vector<ScalarT, IdxT>* y,
                                                         Vector<ScalarT, IdxT>* x)
    {
      // x = beta*x + alpha*V*y OR x = beta*x + alpha*V^Ty
      const ScalarT* V_data = V->getData(memory::HOST);
      const ScalarT* y_data = y->getData(memory::HOST);
      ScalarT*       x_data = x->getData(memory::HOST);
      const IdxT     n      = V->getSize();

      IdxT    i, j;
      ScalarT sum;
      switch (transpose)
      {
      case 'T':
        assert((V->getSize() == y->getSize())
               && "gemv: Size mismatch! Size of V does not match size of y.");
        for (i = 0; i < k; ++i)
        {
          sum       = beta * x_data[i];
          ScalarT c = 0.0;
          for (j = 0; j < n; ++j)
          {
            ScalarT y = (alpha * V_data[i * n + j] * y_data[j]) - c;
            ScalarT t = sum + y;
            c         = (t - sum) - y;
            sum       = t;
          }
          x_data[i] = sum;
        }
        break;
      case 'N':
        assert((V->getSize() == x->getSize())
               && "gemv: Size mismatch! Size of V does not match size of x.");
        for (i = 0; i < n; ++i)
        {
          sum       = beta * x_data[i];
          ScalarT c = 0.0;
          for (j = 0; j < k; ++j)
          {
            ScalarT y = (alpha * V_data[n * j + i] * y_data[j]) - c;
            ScalarT t = sum + y;
            c         = (t - sum) - y;
            sum       = t;
          }
          x_data[i] = sum;
        }
        break;
      default:
        out::error() << "Unrecognized transpose option " << transpose
                     << " in gemv. Valid options are 'N' (not transposed) and 'T' (transposed).\n";
      } // switch
      x->setDataUpdated(memory::HOST);
    }

    /**
     * @brief mass (bulk) axpy i.e, y = y - x*alpha where alpha is a vector
     *
     * @param[in] size number of elements in y
     * @param[in] alpha vector size k x 1
     * @param[in] x (multi)vector size size x k
     * @param[in,out] y vector size size x 1 (this is where the result is stored)
     *
     * @pre _k_ > 0, _size_ > 0, _size_ = x->getSize()
     */
    template <typename scalar_type, typename index_type>
    void VectorHandlerCpu<scalar_type, index_type>::axpyMulti(IdxT                   size,
                                                              Vector<ScalarT, IdxT>* alpha,
                                                              IdxT                   k,
                                                              Vector<ScalarT, IdxT>* x,
                                                              Vector<ScalarT, IdxT>* y)
    {
      ScalarT* alpha_data = alpha->getData(memory::HOST);
      ScalarT* y_data     = y->getData(memory::HOST);
      ScalarT* x_data     = x->getData(memory::HOST);
      IdxT     i, j;
      ScalarT  sum;

      for (i = 0; i < size; ++i)
      {
        sum = 0.0;
        for (j = 0; j < k; ++j)
        {
          sum += x_data[j * size + i] * alpha_data[j];
        }
        y_data[i] = y_data[i] - sum;
      }
      y->setDataUpdated(memory::HOST);
    }

    /**
     * @brief mass (bulk) dot product i.e, V^T x, where V is n x k dense multivector (a dense
     * multivector consisting of k vectors size n) and x is k x 2 dense multivector (a multivector
     * consisting of two vectors size n each)
     *
     * @param[in] size Number of elements in a single vector in V
     * @param[in] V Multivector; k vectors size n x 1 each
     * @param[in] k Number of vectors in V
     * @param[in] x Multivector; 2 vectors size n x 1 each
     * @param[out] res Multivector; 2 vectors size k x 1 each (result is returned in res)
     *
     * @pre _size_ > 0, _k_ > 0, size = x->getSize(), _res_ needs to be allocated
     */
    template <typename scalar_type, typename index_type>
    void VectorHandlerCpu<scalar_type, index_type>::dot2Multi(IdxT                   size,
                                                              Vector<ScalarT, IdxT>* V,
                                                              IdxT                   k,
                                                              Vector<ScalarT, IdxT>* x,
                                                              Vector<ScalarT, IdxT>* res)
    {
      ScalarT*       res_data = res->getData(memory::HOST);
      const ScalarT* x_data   = x->getData(memory::HOST);
      const ScalarT* V_data   = V->getData(memory::HOST);

      ScalarT c0 = 0.0;
      ScalarT cq = 0.0;

      for (IdxT i = 0; i < k; ++i)
      {
        res_data[i]     = 0.0;
        res_data[i + k] = 0.0;

        // Make sure we don't accumulate round-off errors
        for (IdxT j = 0; j < size; ++j)
        {
          ScalarT y0 = (V_data[i * size + j] * x_data[j]) - c0;
          ScalarT yq = (V_data[i * size + j] * x_data[j + size]) - cq;
          ScalarT t0 = res_data[i] + y0;
          ScalarT tq = res_data[i + k] + yq;
          c0         = (t0 - res_data[i]) - y0;
          cq         = (tq - res_data[i + k]) - yq;

          res_data[i]     = t0;
          res_data[i + k] = tq;
        }
      }
      res->setDataUpdated(memory::HOST);
    }

    /**
     * @brief Scale a vector by a diagonal matrix
     *
     * @param[in] diag Diagonal vector
     * @param[in,out] vec Vector to be scaled
     */
    template <typename scalar_type, typename index_type>
    void VectorHandlerCpu<scalar_type, index_type>::scal(Vector<ScalarT, IdxT>* diag, Vector<ScalarT, IdxT>* vec)
    {
      const ScalarT* diag_data = diag->getData(memory::HOST);
      ScalarT*       vec_data  = vec->getData(memory::HOST);
      IdxT           n         = vec->getSize();

      for (IdxT i = 0; i < n; ++i)
      {
        vec_data[i] *= diag_data[i];
      }
      vec->setDataUpdated(memory::HOST);
    }

    /**
     * @brief Scale a vector by a diagonal matrix represented by a contiguous subvector of diag
     *
     * @param[in] diag Diagonal vector
     * @param[in,out] vec Vector to be scaled
     * @param[in] diag_offset - the index of diag where the diagonal matrix begins
     */
    template <typename scalar_type, typename index_type>
    void VectorHandlerCpu<scalar_type, index_type>::scal(Vector<ScalarT, IdxT>* diag, Vector<ScalarT, IdxT>* vec, IdxT diag_offset)
    {
      const ScalarT* diag_data = &diag->getData(memory::HOST)[diag_offset];
      ScalarT*       vec_data  = vec->getData(memory::HOST);
      IdxT           n         = vec->getSize();

      for (IdxT i = 0; i < n; ++i)
      {
        vec_data[i] *= diag_data[i];
      }
      vec->setDataUpdated(memory::HOST);
    }

    /**
     * @brief Multiplies vector by an inverse of a diagonal matrix.
     *  This is equivalent to solving a system with the original matrix,
     *  with the vector as the right hand side (typically used in this context).
     *
     *
     * @param[in]  diag   - diagonal matrix stored in a vector object
     * @param[in,out] vec - vector to be divided
     *
     * @pre The two vectors must be the same size
     *
     * @return 0 if successful, 1 otherwise
     */
    template <typename scalar_type, typename index_type>
    int VectorHandlerCpu<scalar_type, index_type>::diagSolve(Vector<ScalarT, IdxT>* diag, Vector<ScalarT, IdxT>* vec)
    {
      ScalarT* diag_data = diag->getData(memory::HOST);
      ScalarT* vec_data  = vec->getData(memory::HOST);
      IdxT     n         = vec->getSize();

      for (IdxT i = 0; i < n; ++i)
      {
        vec_data[i] /= diag_data[i];
      }
      vec->setDataUpdated(memory::HOST);
      return 0;
    }

    /**
     * @brief Take the element-wise max of two vectors.
     * Each element of the output will be the maximum value of the corresponding elements in the
     * input vectors.
     *
     * @param[in]  x   - First input vector
     * @param[in]  y   - Second input vector
     * @param[out] out - Output vector
     *
     * @pre The three vectors must be the same size
     *
     * @return 0 if successful, 1 otherwise
     */
    template <typename scalar_type, typename index_type>
    int VectorHandlerCpu<scalar_type, index_type>::max(Vector<ScalarT, IdxT>* x, Vector<ScalarT, IdxT>* y, Vector<ScalarT, IdxT>* out)
    {
      const ScalarT* x_data   = x->getData(memory::HOST);
      const ScalarT* y_data   = y->getData(memory::HOST);
      ScalarT*       out_data = out->getData(memory::HOST);
      IdxT           n        = y->getSize();

      for (IdxT i = 0; i < n; ++i)
      {
        out_data[i] = std::max(x_data[i], y_data[i]);
      }
      out->setDataUpdated(memory::HOST);
      return 0;
    }

    /**
     * @brief Take the element-wise absolute value of a vector.
     *
     * @param[in]  in  - Input vector
     * @param[out] out - Output vector
     *
     * @return 0 if successful, 1 otherwise
     */
    template <typename scalar_type, typename index_type>
    int VectorHandlerCpu<scalar_type, index_type>::abs(Vector<ScalarT, IdxT>* in, Vector<ScalarT, IdxT>* out)
    {
      const ScalarT* in_data  = in->getData(memory::HOST);
      ScalarT*       out_data = out->getData(memory::HOST);
      IdxT           n        = in->getSize();

      for (IdxT i = 0; i < n; ++i)
      {
        out_data[i] = std::abs(in_data[i]);
      }
      out->setDataUpdated(memory::HOST);
      return 0;
    }

    // Available template instantiations
    template class VectorHandlerCpu<double, size_t>;
    template class VectorHandlerCpu<double, int>;

  } // namespace LinearAlgebra
} // namespace GridKit
