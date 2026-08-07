/**
 * @file EnzymeTests.hpp
 * @author Nicholson Koukpaizan (koukpaizannk@ornl.gov)
 *
 */

#include <iomanip>
#include <iostream>

#include <GridKit/AutomaticDifferentiation/Enzyme/EnzymeDefinitions.hpp>
#include <GridKit/AutomaticDifferentiation/Enzyme/LowerSparseStorage.hpp>
#include <GridKit/LinearAlgebra/SparseMatrix/CooMatrix.hpp>
#include <GridKit/MemoryUtilities/MemoryUtils.hpp>
#include <GridKit/Testing/TestHelpers.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    template <class ScalarT, typename IdxT>
    class EnzymeTests
    {
    public:
      using SparseMatrix = GridKit::LinearAlgebra::CooMatrix<ScalarT, IdxT>;

      EnzymeTests()  = default;
      ~EnzymeTests() = default;

      /**
       * @brief Standalone Enzyme test that computes the derivative of a scalar function (square)
       *
       * The function to differentiate is f(x) = x^2
       * The expected derivative if f'(x) = 2*x
       *
       */
      TestOutcome scalarFunction()
      {
        TestStatus success = true;

        using namespace GridKit::Enzyme::Scalar;

        ScalarT var        = 5.0;
        ScalarT sq         = square(var);
        ScalarT dsq_ref    = dsquare(var);
        ScalarT dsq_enzyme = __enzyme_autodiff<ScalarT>(square, var);

        success *= GridKit::Testing::isEqual(dsq_enzyme, dsq_ref);

        if (!success)
        {
          std::cout << "x = " << var << "\n"
                    << "x^2 = " << sq << "\n"
                    << "Reference d(x^2)/dx = " << dsq_ref << "\n"
                    << "Enzyme d(x^2)/dx = " << dsq_enzyme << "\n";
        }

        return success.report(__func__);
      }

      /**
       * @brief Standalone Enzyme test that computes the derivative of a vector function
       *
       * The function to differentiate is f_i(x) = sum(x_j^2) for j <= i
       * The expected Jacobian is a lower triangular matrix of f_{i,j}'(x) = 2*x_j if j <= i
       *
       */
      TestOutcome vectorFunctionSparse()
      {
        TestStatus success = true;

        IdxT     N = 5;
        ScalarT* x = new ScalarT[static_cast<size_t>(N)];
        ScalarT* y = new ScalarT[static_cast<size_t>(N)];

        // Input initialization to [0.0, 1.0, 2.0, ...]
        for (IdxT i = 0; i < N; ++i)
        {
          x[i] = static_cast<ScalarT>(i);
        }

        // Output function evaluation
        vectorFunction(N, x, y);

        // Enzyme Jacobian
        SparseMatrix* jac_enzyme = vectorFunctionSparseEnzymeJacobian(N, x);

        // Reference Jacobian
        SparseMatrix* jac_ref = vectorFunctionSparseReferenceJacobian(N, x);

        // Compare the Jacobians
        success *= checkVectorFunctionSparseJacobian(jac_ref, jac_enzyme);

        if (!success)
        {
          std::cout << "Enzyme Jacobian\n";
          jac_enzyme->print();
          std::cout << "Reference Jacobian\n";
          jac_ref->print();
        }

        delete[] x;
        delete[] y;
        delete jac_enzyme;
        delete jac_ref;

        return success.report(__func__);
      }

    private:
      static ScalarT square(ScalarT x)
      {
        return x * x;
      }

      static ScalarT dsquare(ScalarT x)
      {
        return 2.0 * x;
      }

      FORCE_INLINE static void vectorFunction(IdxT N, ScalarT* x, ScalarT* y)
      {
        for (IdxT idx = 0; idx < N; ++idx)
        {
          y[idx] = 0.0;
          for (IdxT idy = 0; idy <= idx; idy++)
          {
            y[idx] += square(x[idy]);
          }
        }
      }

      SparseMatrix* vectorFunctionSparseEnzymeJacobian(const IdxT N, const ScalarT* x)
      {
        using namespace GridKit::Enzyme::Sparse;

        IdxT* index_maps = new IdxT[static_cast<size_t>(N)];
        for (IdxT i = 0; i < N; i++)
        {
          index_maps[i] = i;
        }

        size_t*  rows_buffer = new size_t[N * N];
        size_t*  cols_buffer = new size_t[N * N];
        ScalarT* vals_buffer = new ScalarT[N * N];
        size_t   current_nnz = 0;
        for (size_t i = 0; i < N; i++)
        {
          ScalarT* output   = __enzyme_todense<ScalarT*>((void*) ident_load<ScalarT, size_t>,
                                                       (void*) ident_store<ScalarT, size_t>,
                                                       i);
          ScalarT* d_output = __enzyme_todense<ScalarT*>((void*) sparse_load<ScalarT, size_t>,
                                                         (void*) sparse_store<ScalarT, size_t>,
                                                         i,
                                                         1.0,
                                                         index_maps,
                                                         index_maps,
                                                         rows_buffer,
                                                         cols_buffer,
                                                         vals_buffer,
                                                         &current_nnz);

          __enzyme_fwddiff<void>((void*) vectorFunction,
                                 enzyme_const,
                                 N,
                                 enzyme_dup,
                                 x,
                                 output,
                                 enzyme_dupnoneed,
                                 (ScalarT*) 0x1,
                                 d_output);
        }

        size_t*  rows = new size_t[current_nnz];
        size_t*  cols = new size_t[current_nnz];
        ScalarT* vals = new ScalarT[current_nnz];
        for (size_t ind = 0; ind < current_nnz; ++ind)
        {
          rows[ind] = rows_buffer[ind];
          cols[ind] = cols_buffer[ind];
          vals[ind] = vals_buffer[ind];
        }

        delete[] index_maps;
        delete[] rows_buffer;
        delete[] cols_buffer;
        delete[] vals_buffer;

        // Hijacking constructor sets rows, cols and vals to nullptr
        return new SparseMatrix(N, N, current_nnz, &rows, &cols, &vals);
      }

      SparseMatrix* vectorFunctionSparseReferenceJacobian(const IdxT N, const ScalarT* x)
      {
        IdxT     nnz         = N * (N + 1) / 2; // lower triangular matrix
        IdxT*    rows        = new IdxT[static_cast<size_t>(nnz)];
        IdxT*    cols        = new IdxT[static_cast<size_t>(nnz)];
        ScalarT* vals        = new ScalarT[static_cast<size_t>(nnz)];
        IdxT     current_nnz = 0;
        for (IdxT idy = 0; idy < N; ++idy)
        {
          for (IdxT idx = 0; idx < N; ++idx)
          {
            if (idy <= idx && current_nnz < nnz)
            {
              rows[current_nnz] = idx;
              cols[current_nnz] = idy;
              vals[current_nnz] = dsquare(x[idy]);
              current_nnz++;
            }
          }
        }

        // Hijacking constructor sets rows, cols and vals to nullptr
        return new SparseMatrix(N, N, current_nnz, &rows, &cols, &vals);
      }

      // @todo move to a common location if this needs to be reused
      TestStatus checkVectorFunctionSparseJacobian(SparseMatrix* matrix_1,
                                                   SparseMatrix* matrix_2)
      {
        TestStatus success = true;

        IdxT     nnz_1  = matrix_1->getNnz();
        IdxT*    rows_1 = matrix_1->getRowData();
        IdxT*    cols_1 = matrix_1->getColData();
        ScalarT* vals_1 = matrix_1->getValues();

        IdxT     nnz_2  = matrix_2->getNnz();
        IdxT*    rows_2 = matrix_2->getRowData();
        IdxT*    cols_2 = matrix_2->getColData();
        ScalarT* vals_2 = matrix_2->getValues();

        success *= (nnz_1 == nnz_2);

        for (IdxT ind = 0; ind < std::min(nnz_1, nnz_2); ++ind)
        {
          success *= (rows_1[ind] == rows_2[ind]);
          success *= (cols_1[ind] == cols_2[ind]);
          success *= GridKit::Testing::isEqual(vals_1[ind], vals_2[ind]);
        }

        return success;
      }
    };

  } // namespace Testing
} // namespace GridKit
