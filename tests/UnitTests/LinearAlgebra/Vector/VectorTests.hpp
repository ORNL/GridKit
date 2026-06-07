#pragma once
#include <algorithm>
#include <iomanip>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

#include <GridKit/LinearAlgebra/MemoryUtils.hpp>
#include <GridKit/LinearAlgebra/Vector/Vector.hpp>
#include <GridKit/Testing/TestHelpers.hpp>

namespace GridKit
{
  namespace Testing
  {
    using namespace LinearAlgebra;

    /**
     * @class Tests for vector operations.
     */
    template <class ScalarT, typename IdxT>
    class VectorTests
    {
    public:
      VectorTests(memory::MemorySpace memspace = memory::HOST)
        : memspace_(memspace)
      {
      }

      virtual ~VectorTests()
      {
      }

      /**
       * @brief Test vector constructor with specified size and number of vectors.
       *
       * @param[in] N Number of elements in the vector.
       * @param[in] k Number of vectors in the multivector.
       * @return TestOutcome indicating success or failure of the test.
       */
      TestOutcome vectorConstructor(IdxT N, IdxT k)
      {
        TestStatus status;
        status = true;

        Vector<ScalarT, IdxT> x(N, k);

        if (x.getCapacity() != N)
        {
          std::cout << "The capacity of the vector is " << x.getCapacity()
                    << ", expected: " << N << "\n";
          status *= false;
        }

        if (x.getSize() != N)
        {
          std::cout << "The size of the vector is " << x.getSize()
                    << ", expected: " << N << "\n";
          status *= false;
        }

        if (x.getNumVectors() != k)
        {
          std::cout << "The number of vectors in the multivector is " << x.getNumVectors()
                    << ", expected: " << k << "\n";
          status *= false;
        }

        return status.report(__func__);
      }

      /**
       * @brief Test vector constructor with specified size and default number of vectors (1).
       *
       * @param[in] N Number of elements in the vector.
       * @return TestOutcome indicating success or failure of the test.
       */
      TestOutcome vectorConstructor(IdxT N)
      {
        return vectorConstructor(N, 1);
      }

      /**
       * @brief Test resizing a vector to a new size.
       *
       * @param[in] N Current size of the vector.
       * @param[in] new_N New size to which the vector should be resized.
       * @return TestOutcome indicating success or failure of the test.
       */
      TestOutcome resize(IdxT N, IdxT new_N)
      {
        TestStatus status;
        status = true;

        Vector<ScalarT, IdxT> x(N);
        x.allocate(memspace_);

        x.resize(new_N);

        if (x.getSize() != new_N)
        {
          std::cout << "The size of the vector after resizing is " << x.getSize()
                    << ", expected: " << new_N << "\n";
          status *= false;
        }

        return status.report(__func__);
      }

      /**
       * @brief Test setting data in a vector from array.
       *
       * @param[in] N Number of elements in the vector.
       * @return TestOutcome indicating success or failure of the test.
       */
      TestOutcome setData(IdxT N)
      {
        TestStatus status;
        status = true;

        Vector<ScalarT, IdxT> x(N);

        ScalarT* data = new ScalarT[N];
        for (IdxT i = 0; i < N; ++i)
        {
          data[i] = 0.1 * (ScalarT) i;
        }
        x.setData(data, memspace_);

        const ScalarT* x_data = x.getData(memspace_);

        if (x_data == nullptr)
        {
          std::cout << "The data pointer is null after setting.\n";
          status *= false;
        }
        else
        {
          for (IdxT i = 0; i < N; ++i)
          {
            if (!isEqual(x_data[i], data[i]))
            {
              std::cout << "The data in the vector is incorrect at index " << i
                        << ", expected: " << data[i]
                        << ", got: " << x_data[i] << "\n";
              status *= false;
              break;
            }
          }
        }

        delete[] data;
        return status.report(__func__);
      }

      /**
       * @brief Test copying data between vector-array and vector-vector.
       *
       * This creates an array, copies it to a vector in the current memory space, then
       * copies it to another vector in the same memory space, and finally back to a third on the
       * HOST. Then, it verifies the content of the final vector. This test only passes if all
       * copies are successful.
       *
       * @param[in] N Number of elements in the vector.
       * @return TestOutcome indicating success or failure of the test.
       */
      TestOutcome copyFromExternal(IdxT N)
      {
        TestStatus status;
        status = true;

        Vector<ScalarT, IdxT> x(N);
        ScalarT*              data = new ScalarT[N];
        for (IdxT i = 0; i < N; ++i)
        {
          data[i] = 0.1 * (ScalarT) i;
        }

        // array -> memspace
        x.allocate(memspace_);
        x.copyFromExternal(data, memory::HOST, memspace_);

        // memspace -> memspace
        Vector<ScalarT, IdxT> y(N);
        y.allocate(memspace_);
        y.copyFromExternal(x, memspace_, memspace_);

        // memspace -> host
        Vector<ScalarT, IdxT> z(N);
        z.allocate(memory::HOST);
        z.copyFromExternal(y, memspace_, memory::HOST);

        const ScalarT* z_data = z.getData(memory::HOST);

        if (z_data == nullptr)
        {
          std::cout << "The data pointer is null after copying from vector.\n";
          status *= false;
        }
        else
        {
          for (IdxT i = 0; i < N; ++i)
          {
            if (!isEqual(z_data[i], data[i]))
            {
              std::cout << "The data in the copied vector is incorrect at index " << i
                        << ", expected: " << data[i]
                        << ", got: " << z_data[i] << "\n";
              status *= false;
              break;
            }
          }
        }

        delete[] data;
        return status.report(__func__);
      }

      /**
       * @brief Test copying data from vector to an array.
       *
       * This creates a vector, copies data to it, and then copies the data to an array
       * in the current memory space. Finally, it uses the MemoryHandler to copy the data
       * to HOST for verification.
       *
       * @param[in] N Number of elements in the vector.
       * @return TestOutcome indicating success or failure of the test.
       */
      TestOutcome copyToExternal(IdxT N)
      {
        TestStatus status = true;

        Vector   x(N);
        ScalarT* data = new ScalarT[N];
        for (int i = 0; i < N; ++i)
        {
          data[i] = 0.1 * (ScalarT) i;
        }

        x.allocate(memspace_);
        x.copyFromExternal(data, memory::HOST, memspace_);

        // Copy data to an array on current memspace
        ScalarT* dest = new ScalarT[N];
        // second argument is source, third is destination
        x.copyToExternal(dest, memspace_, memspace_);

        // Copy to host to verify
        ScalarT* dest_h = new ScalarT[N];
        if (memspace_ == memory::DEVICE)
        {
          mh_.copyArrayDeviceToHost(dest_h, dest, N);
          dest = dest_h;
        }
        else
        {
          // If we are on HOST, we can use dest directly
          delete[] dest_h;
          dest_h = dest;
        }

        // Verify the copied data
        for (int i = 0; i < N; ++i)
        {
          if (!isEqual(dest_h[i], data[i]))
          {
            std::cout << "The data in the destination array is incorrect at index " << i
                      << ", expected: " << data[i]
                      << ", got: " << dest_h[i] << "\n";
            status *= false;
            break;
          }
        }

        delete[] data;
        delete[] dest;
        return status.report(__func__);
      }

      /**
       * @brief Test setting all elements of a vector to a constant value.
       *
       * @param[in] N Number of elements in the vector.
       * @return TestOutcome indicating success or failure of the test.
       */
      TestOutcome setToConst(IdxT N)
      {
        constexpr ScalarT ONE  = 1.0;
        constexpr ScalarT ZERO = 0.0;

        TestStatus success = true;

        IdxT vector_size    = N;
        IdxT number_vectors = 3;

        Vector<ScalarT, IdxT> x(vector_size, number_vectors);
        x.allocate(memspace_);
        if (memspace_ == memory::DEVICE)
          x.allocate(memory::HOST);

        x.setToZero(memspace_);
        success *= verifyAnswer(x, ZERO);

        x.setToConst(1, ONE, memspace_); // set vector 1 to ones
        success *= verifyAnswer(vector_size, x.getData(0, memspace_), ZERO);
        success *= verifyAnswer(vector_size, x.getData(1, memspace_), ONE);
        success *= verifyAnswer(vector_size, x.getData(2, memspace_), ZERO);

        x.setToConst(ONE, memspace_);
        success *= verifyAnswer(x, ONE);

        x.setToZero(1, memspace_); // set vector 1 to zeros
        success *= verifyAnswer(vector_size, x.getData(0, memspace_), ONE);
        success *= verifyAnswer(vector_size, x.getData(1, memspace_), ZERO);
        success *= verifyAnswer(vector_size, x.getData(2, memspace_), ONE);

        return success.report(__func__);
      }

      /**
       * @brief Test syncing data between HOST and DEVICE memory spaces.
       *
       * Creates a vector allocated in the specified memory space, then sync
       * to the other memory space and verify that the data is synced correctly.
       *
       * @param[in] N Number of elements in the vector.
       * @return TestOutcome returns a report on the test.
       */
      TestOutcome syncData(IdxT N = 4)
      {
        constexpr ScalarT ONE  = 1.0;
        constexpr ScalarT ZERO = 0.0;

        TestStatus success;
        success = true;

        if (memspace_ == memory::HOST)
        {
          return success.report(__func__);
        }

        IdxT vector_size    = N;
        IdxT number_vectors = 3;

        Vector x(vector_size, number_vectors);
        x.allocate(memspace_);
        x.allocate(memory::HOST);

        // Set all vectors in x on device to ones
        x.setToConst(ONE, memspace_);
        // Sync host (all ones on the host, as well)
        x.syncData(memory::HOST);
        // Set vector 1 to all zeros on host
        x.setToZero(1, memory::HOST);
        // Sync vector 1 on device
        x.syncData(1, memspace_);

        // Check what we have on device now is correct
        success *= verifyAnswer(vector_size, x.getData(0, memspace_), ONE);
        success *= verifyAnswer(vector_size, x.getData(1, memspace_), ZERO);
        success *= verifyAnswer(vector_size, x.getData(2, memspace_), ONE);

        return success.report(__func__);
      }

    private:
      memory::MemorySpace memspace_{memory::HOST};
      MemoryHandler       mh_;

      /// Check if vector elements are set to the same number
      bool verifyAnswer(Vector<ScalarT, IdxT>& x, ScalarT answer)
      {
        bool success = true;

        if (memspace_ == memory::DEVICE)
        {
          x.syncData(memory::HOST);
        }

        for (IdxT i = 0; i < x.getSize(); ++i)
        {
          // std::cout << x->getData("cpu")[i] << "\n";
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

      /// Check if an array elements are set to the same number
      bool verifyAnswer(IdxT size, const ScalarT* data, ScalarT answer)
      {
        bool     success = true;
        ScalarT* x       = nullptr;

        // If the data is on device copy it to the host
        if (memspace_ == memory::DEVICE)
        {
          mh_.allocateArrayOnHost(&x, size);
          mh_.copyArrayDeviceToHost(x, data, size);
          // Set `data` to point to the host copy
          data = x;
        }

        for (size_t i = 0; i < static_cast<size_t>(size); ++i)
        {
          if (!isEqual(data[i], answer))
          {
            std::cout << std::setprecision(16);
            success = false;
            std::cout << "Solution vector element x[" << i << "] = " << data[i]
                      << ", expected: " << answer << "\n";
            break;
          }
        }

        if (memspace_ == memory::DEVICE)
        {
          mh_.deleteOnHost(x);
          data = nullptr;
        }

        return success;
      }
    }; // class
  } // namespace Testing
} // namespace GridKit
