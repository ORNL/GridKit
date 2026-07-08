#pragma once
#include <cmath>
#include <iomanip>
#include <memory>

#include <GridKit/LinearAlgebra/Vector/Vector.hpp>
#include <GridKit/LinearAlgebra/Vector/VectorHandler.hpp>
#include <GridKit/MemoryUtilities/MemoryUtils.hpp>
#include <GridKit/Testing/TestHelpers.hpp>

namespace GridKit
{
  namespace Testing
  {
    using namespace LinearAlgebra;

    /**
     * @class Tests for the vector handler.
     */
    template <class ScalarT, typename IdxT>
    class VectorHandlerTests
    {
    public:
      VectorHandlerTests(VectorHandler<ScalarT, IdxT>& handler, memory::MemorySpace memspace = memory::HOST)
        : handler_(handler),
          memspace_(memspace)
      {
      }

      virtual ~VectorHandlerTests()
      {
      }

      TestOutcome vectorHandlerConstructor()
      {
        TestStatus status;
        status.skipTest();

        return status.report(__func__);
      }

      /**
       * @brief Test vector infinity norm.
       */
      TestOutcome amax(IdxT N)
      {
        TestStatus status;
        status = true;

        Vector<ScalarT, IdxT> x(N);

        ScalarT* data = new ScalarT[N];
        for (IdxT i = 0; i < N; ++i)
        {
          data[i] = static_cast<ScalarT>(0.1) * static_cast<ScalarT>(i);
        }
        x.allocate(memspace_);
        x.copyFromExternal(data, memory::HOST, memspace_);

        ScalarT result = handler_.amax(&x, memspace_);
        ScalarT answer = static_cast<ScalarT>(N - 1) * static_cast<ScalarT>(0.1);

        if (!isEqual(result, answer))
        {
          std::cout << "The result " << result << " is incorrect. "
                    << "Expected answer is " << answer << "\n";
          status *= false;
        }

        delete[] data;
        return status.report(__func__);
      }

      /**
       * @brief Test axpy: y = alpha*x + y.
       */
      TestOutcome axpy(IdxT N)
      {
        TestStatus status;
        status = true;

        Vector<ScalarT, IdxT> x(N);
        Vector<ScalarT, IdxT> y(N);

        x.allocate(memspace_);
        y.allocate(memspace_);

        x.setToConst(3.0, memspace_);
        y.setToConst(1.0, memspace_);

        ScalarT alpha = static_cast<ScalarT>(0.5);

        // the result is a vector with y[i] = 2.5 for all i;
        handler_.axpy(alpha, &x, &y, memspace_);
        status *= verifyAnswer(y, 2.5);

        return status.report(__func__);
      }

      /**
       * @brief Test dot product.
       */
      TestOutcome dot(IdxT N)
      {
        TestStatus status;
        status = true;

        Vector<ScalarT, IdxT> x(N);
        Vector<ScalarT, IdxT> y(N);

        x.allocate(memspace_);
        y.allocate(memspace_);

        x.setToConst(0.25, memspace_);
        y.setToConst(4.0, memspace_);

        // the answer is N
        ScalarT answer = static_cast<ScalarT>(N);
        ScalarT result = handler_.dot(&x, &y, memspace_);

        if (!isEqual(result, answer))
        {
          std::cout << "The result " << result << " is incorrect. "
                    << "Expected answer is " << answer << "\n";
          status *= false;
        }

        return status.report(__func__);
      }

      /**
       * @brief Test scal: x = alpha*x.
       */
      TestOutcome scal(IdxT N)
      {
        TestStatus status;
        status = true;

        Vector<ScalarT, IdxT> x(N);

        x.allocate(memspace_);

        x.setToConst(1.25, memspace_);

        ScalarT alpha = 3.5;

        // the answer is x[i] = 4.375;
        ScalarT answer = 4.375;
        handler_.scal(alpha, &x, memspace_);
        status *= verifyAnswer(x, answer);

        return status.report(__func__);
      }

      /**
       * @brief Test mass (bulk) axpy: y = y - x*alpha, where x is [N x K] and alpha is [K x 1].
       */
      TestOutcome axpyMulti(IdxT N, IdxT K)
      {
        TestStatus status;
        status = true;

        Vector<ScalarT, IdxT> x(N, K);
        Vector<ScalarT, IdxT> y(N);
        Vector<ScalarT, IdxT> alpha(K);

        x.allocate(memspace_);
        y.allocate(memspace_);
        alpha.allocate(memspace_);

        alpha.setToConst(-1.0, memspace_);
        y.setToConst(2.0, memspace_);

        for (IdxT ii = 0; ii < K; ++ii)
        {
          ScalarT c;
          if (ii % 2 == 0)
          {
            c = -1.0;
          }
          else
          {
            c = static_cast<ScalarT>(0.5);
          }
          x.setToConst(ii, c, memspace_);
        }

        IdxT    r   = K % 2;
        ScalarT res = static_cast<ScalarT>((std::floor(static_cast<double>(K) / 2.0) + static_cast<double>(r)) * 1.0
                                           + std::floor(static_cast<double>(K) / 2.0) * (-0.5));

        handler_.axpyMulti(N, &alpha, K, &x, &y, memspace_);
        status *= verifyAnswer(y, 2.0 - res);

        return status.report(__func__);
      }

      /**
       * @brief Test mass (bulk) dot product: V^T x.
       */
      TestOutcome massDot(IdxT N, IdxT K)
      {
        TestStatus status;
        status = true;

        Vector<ScalarT, IdxT> x(N, K);
        Vector<ScalarT, IdxT> y(N, 2);
        Vector<ScalarT, IdxT> res(K, 2);
        x.allocate(memspace_);
        y.allocate(memspace_);
        res.allocate(memspace_);
        // res is a write-only output; mark it up to date before it is used.
        res.setToZero(memspace_);

        x.setToConst(1.0, memspace_);
        y.setToConst(-1.0, memspace_);
        handler_.dot2Multi(N, &x, K, &y, &res, memspace_);

        status *= verifyAnswer(res, (-1.0) * static_cast<ScalarT>(N));

        return status.report(__func__);
      }

      /**
       * @brief Test dense matrix-vector product (gemv), transposed and not.
       */
      TestOutcome gemv(IdxT N, IdxT K)
      {
        TestStatus status;
        status = true;

        Vector<ScalarT, IdxT> V(N, K);
        Vector<ScalarT, IdxT> yN(K); ///< For the test with NO TRANSPOSE
        Vector<ScalarT, IdxT> xN(N);
        Vector<ScalarT, IdxT> yT(N); ///< for the test with TRANSPOSE
        Vector<ScalarT, IdxT> xT(K);

        V.allocate(memspace_);
        yN.allocate(memspace_);
        xN.allocate(memspace_);
        yT.allocate(memspace_);
        xT.allocate(memspace_);

        V.setToConst(1.0, memspace_);
        yN.setToConst(-1.0, memspace_);
        xN.setToConst(0.5, memspace_);
        yT.setToConst(-1.0, memspace_);
        xT.setToConst(0.5, memspace_);

        ScalarT alpha = -1.0;
        ScalarT beta  = 1.0;
        handler_.gemv('N', K, alpha, beta, &V, &yN, &xN, memspace_);
        status *= verifyAnswer(xN, static_cast<ScalarT>(K) + static_cast<ScalarT>(0.5));
        handler_.gemv('T', K, alpha, beta, &V, &yT, &xT, memspace_);
        status *= verifyAnswer(xT, static_cast<ScalarT>(N) + static_cast<ScalarT>(0.5));

        return status.report(__func__);
      }

      /**
       * @brief Test scaling a vector by a diagonal matrix.
       */
      TestOutcome scale(IdxT N)
      {
        TestStatus status;
        status = true;

        Vector<ScalarT, IdxT> diag(N);
        Vector<ScalarT, IdxT> vec(N);

        // diag[i] = i + 1, vec[i] = 3.0
        // expected result vec[i] = (i + 1) * 3.0
        diag.allocate(memspace_);
        vec.allocate(memspace_);

        vec.setToConst(3.0, memspace_);

        auto diag_data = std::unique_ptr<ScalarT[]>(new ScalarT[N]);
        for (IdxT i = 0; i < N; ++i)
        {
          diag_data[i] = static_cast<ScalarT>(i + 1);
        }
        diag.copyFromExternal(diag_data.get(), memory::HOST, memspace_);

        handler_.scal(&diag, &vec, memspace_);

        if (memspace_ == memory::DEVICE)
        {
          vec.syncData(memory::HOST);
        }

        for (IdxT i = 0; i < N; ++i)
        {
          if (!isEqual(vec.getData(memory::HOST)[i], static_cast<ScalarT>(i + 1) * 3.0))
          {
            std::cout << "Solution vector element vec[" << i << "] = " << vec.getData(memory::HOST)[i]
                      << ", expected: " << static_cast<ScalarT>(i + 1) * 3.0 << "\n";
            status *= false;
            break;
          }
        }

        return status.report(__func__);
      }

      /**
       * @brief Test dividing a vector by a diagonal matrix.
       */
      TestOutcome diagSolve(IdxT N)
      {
        TestStatus status;
        status = true;

        Vector<ScalarT, IdxT> diag(N);
        Vector<ScalarT, IdxT> vec(N);

        // diag[i] = i + 1, vec[i] = 3.0
        // expected result vec[i] = 3.0 / (i + 1)
        diag.allocate(memspace_);
        vec.allocate(memspace_);

        vec.setToConst(3.0, memspace_);

        auto diag_data = std::unique_ptr<ScalarT[]>(new ScalarT[N]);
        for (IdxT i = 0; i < N; ++i)
        {
          diag_data[i] = static_cast<ScalarT>(i + 1);
        }
        diag.copyFromExternal(diag_data.get(), memory::HOST, memspace_);

        handler_.diagSolve(&diag, &vec, memspace_);

        if (memspace_ == memory::DEVICE)
        {
          vec.syncData(memory::HOST);
        }

        for (IdxT i = 0; i < N; ++i)
        {
          if (!isEqual(vec.getData(memory::HOST)[i], static_cast<ScalarT>(3.0) / static_cast<ScalarT>(i + 1)))
          {
            std::cout << "Solution vector element vec[" << i << "] = " << vec.getData(memory::HOST)[i]
                      << ", expected: " << static_cast<ScalarT>(3.0) / static_cast<ScalarT>(i + 1) << "\n";
            status *= false;
            break;
          }
        }

        return status.report(__func__);
      }

      /**
       * @brief Test element-wise max of two vectors.
       */
      TestOutcome max(IdxT N)
      {
        TestStatus status;
        status = true;

        Vector<ScalarT, IdxT> x(N);
        Vector<ScalarT, IdxT> y(N);
        Vector<ScalarT, IdxT> z(N);

        x.allocate(memspace_);
        y.allocate(memspace_);
        z.allocate(memspace_);
        // z is a write-only output; mark it up to date before it is used.
        z.setToZero(memspace_);

        auto x_data = std::unique_ptr<ScalarT[]>(new ScalarT[N]);
        auto y_data = std::unique_ptr<ScalarT[]>(new ScalarT[N]);
        for (IdxT i = 0; i < N; ++i)
        {
          if (i % 3 == 0)
          {
            x_data[i] = static_cast<ScalarT>(i + 1);
            y_data[i] = static_cast<ScalarT>(i) * 0.5;
          }
          else
          {
            x_data[i] = -static_cast<ScalarT>(i + 1);
            y_data[i] = static_cast<ScalarT>(i + 1);
          }
        }
        x.copyFromExternal(x_data.get(), memory::HOST, memspace_);
        y.copyFromExternal(y_data.get(), memory::HOST, memspace_);

        handler_.max(&x, &y, &z, memspace_);
        handler_.max(&x, &y, &y, memspace_);

        if (memspace_ == memory::DEVICE)
        {
          y.syncData(memory::HOST);
          z.syncData(memory::HOST);
        }

        for (IdxT i = 0; i < N; ++i)
        {
          if (!isEqual(y.getData(memory::HOST)[i], static_cast<ScalarT>(i + 1)))
          {
            std::cout << "Solution vector element y[" << i << "] = " << y.getData(memory::HOST)[i]
                      << ", expected: " << static_cast<ScalarT>(i + 1) << "\n";
            status *= false;
            break;
          }

          if (!isEqual(z.getData(memory::HOST)[i], static_cast<ScalarT>(i + 1)))
          {
            std::cout << "Solution vector element z[" << i << "] = " << z.getData(memory::HOST)[i]
                      << ", expected: " << static_cast<ScalarT>(i + 1) << "\n";
            status *= false;
            break;
          }
        }

        return status.report(__func__);
      }

      /**
       * @brief Test element-wise absolute value of a vector.
       */
      TestOutcome abs(IdxT N)
      {
        TestStatus status;
        status = true;

        Vector<ScalarT, IdxT> x(N);
        Vector<ScalarT, IdxT> y(N);

        x.allocate(memspace_);
        y.allocate(memspace_);
        // y is a write-only output; mark it up to date before it is used.
        y.setToZero(memspace_);

        auto x_data = std::unique_ptr<ScalarT[]>(new ScalarT[N]);
        for (IdxT i = 0; i < N; ++i)
        {
          if (i % 3 == 0)
          {
            x_data[i] = -static_cast<ScalarT>(i);
          }
          else
          {
            x_data[i] = static_cast<ScalarT>(i);
          }
        }
        x.copyFromExternal(x_data.get(), memory::HOST, memspace_);

        handler_.abs(&x, &y, memspace_);
        handler_.abs(&x, &x, memspace_);

        if (memspace_ == memory::DEVICE)
        {
          x.syncData(memory::HOST);
          y.syncData(memory::HOST);
        }

        for (IdxT i = 0; i < N; ++i)
        {
          if (!isEqual(x.getData(memory::HOST)[i], static_cast<ScalarT>(i)))
          {
            std::cout << "Solution vector element x[" << i << "] = " << x.getData(memory::HOST)[i]
                      << ", expected: " << static_cast<ScalarT>(i) << "\n";
            status *= false;
            break;
          }

          if (!isEqual(y.getData(memory::HOST)[i], static_cast<ScalarT>(i)))
          {
            std::cout << "Solution vector element y[" << i << "] = " << y.getData(memory::HOST)[i]
                      << ", expected: " << static_cast<ScalarT>(i) << "\n";
            status *= false;
            break;
          }
        }

        return status.report(__func__);
      }

    private:
      VectorHandler<ScalarT, IdxT>& handler_;
      memory::MemorySpace           memspace_{memory::HOST};

      // We could verify through norm but that would defeat the purpose of testing vector handler.
      bool verifyAnswer(Vector<ScalarT, IdxT>& x, ScalarT answer)
      {
        bool success = true;

        if (memspace_ == memory::DEVICE)
        {
          x.syncData(memory::HOST);
        }

        for (IdxT i = 0; i < x.getSize(); ++i)
        {
          if (!isEqual(x.getData(memory::HOST)[i], answer))
          {
            std::cout << std::setprecision(16);
            success = false;
            std::cout << "Solution vector element x[" << i << "] = " << x.getData(memory::HOST)[i]
                      << ", expected: " << answer << "\n";
            break;
          }
        }
        return success;
      }
    }; // class
  } // namespace Testing
} // namespace GridKit
