#include <GridKit/LinearAlgebra/SparseMatrix/CooMatrix.hpp>
#include <GridKit/Testing/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    using namespace LinearAlgebra;

    template <class ScalarT, typename IdxT>
    class SparseCooTests
    {

      using CooMatrix = LinearAlgebra::CooMatrix<ScalarT, IdxT>;

    public:
      SparseCooTests(memory::MemorySpace memspace = memory::HOST)
        : memspace_(memspace)
      {
      }

      ~SparseCooTests()
      {
      }

      /**
       * @brief Constructor test for COO sparse matrix
       *
       * Constructs a COO matrix with the given parameters and confirms that the
       * parameters are stored correctly.
       *
       * @param[in] n   - number of rows
       * @param[in] m   - number of columns
       * @param[in] nnz - number of non-zeros
       *
       * @return TestOutcome indicating success or failure of the test
       */
      TestOutcome constructor(IdxT n, IdxT m, IdxT nnz)
      {
        TestStatus success;
        success = true;

        CooMatrix A(n, m, nnz);

        if (A.getNumRows() != n || A.getNumColumns() != m || A.getNnz() != nnz)
        {
          std::cout << "Matrix dimensions do not match expected values.\n";
          success = false;
        }

        return success.report(__func__);
      }

      /**
       * @brief Test set data pointers for the COO sparse matrix
       *
       * Sets the row, column, and value data pointers for the COO matrix and
       * checks if the pointers are set correctly.
       *
       * @param[in] n           - number of rows
       * @param[in] m           - number of columns
       * @param[in] row_density - number of non-zeros per row
       *
       * @return TestOutcome indicating success or failure of the test
       */
      TestOutcome setDataPointers(IdxT n, IdxT m, IdxT row_density)
      {
        TestStatus success;
        success = true;

        IdxT nnz = n * row_density;

        CooMatrix A(n, m, nnz);

        IdxT* h_row_data = new IdxT[static_cast<size_t>(nnz)];
        for (IdxT i = 0; i < nnz; ++i)
        {
          h_row_data[i] = i / row_density;
        }
        IdxT* h_col_data = new IdxT[static_cast<size_t>(nnz)];
        for (IdxT i = 0; i < nnz; ++i)
        {
          h_col_data[i] = i % m;
        }
        ScalarT* h_val_data = new ScalarT[static_cast<size_t>(nnz)];
        for (IdxT i = 0; i < nnz; ++i)
        {
          h_val_data[i] = static_cast<ScalarT>(i + 1);
        }

        if (A.setDataPointers(h_row_data, h_col_data, h_val_data, memory::HOST) != 0)
        {
          std::cout << "Failed to set data pointers.\n";
          success = false;
        }
        else if (A.getRowData(memory::HOST) != h_row_data || A.getColData(memory::HOST) != h_col_data || A.getValues(memory::HOST) != h_val_data)
        {
          std::cout << "Data pointers do not point to expected values.\n";
          success = false;
        }

        delete[] h_row_data;
        delete[] h_col_data;
        delete[] h_val_data;

        return success.report(__func__);
      }

      /**
       * @brief Test getCsrRowData on a small COO matrix with duplicates
       *
       * Builds a 3x3 COO matrix with 5 entries (including one duplicate pair),
       * calls getCsrRowData(), and verifies that the returned CSR row pointers
       * and the updated nnz are correct.
       *
       * @return TestOutcome indicating success or failure of the test
       */
      TestOutcome getCsrRowData()
      {
        TestStatus success;
        success = true;

        const IdxT n       = 3;
        const IdxT m       = 3;
        const IdxT nnz_dup = 5;

        CooMatrix A(n, m, nnz_dup);

        IdxT*    rows = new IdxT[nnz_dup]{2, 0, 0, 1, 1};
        IdxT*    cols = new IdxT[nnz_dup]{2, 0, 1, 1, 1};
        ScalarT* vals = new ScalarT[nnz_dup]{5, 1, 2, 3, 4};

        if (A.setDataPointers(rows, cols, vals, memory::HOST) != 0)
        {
          std::cout << "Failed to set data pointers.\n";
          success = false;
        }
        else
        {
          IdxT* row_ptrs = A.getCsrRowData();

          const IdxT expected_nnz         = 4;
          const IdxT expected_row_ptrs[4] = {0, 2, 3, 4};

          if (A.getNnz() != expected_nnz)
          {
            std::cout << "nnz after deduplication is " << A.getNnz()
                      << ", expected " << expected_nnz << ".\n";
            success = false;
          }

          for (IdxT i = 0; i <= n; ++i)
          {
            if (row_ptrs[i] != expected_row_ptrs[i])
            {
              std::cout << "row_ptrs[" << i << "] = " << row_ptrs[i]
                        << ", expected " << expected_row_ptrs[i] << ".\n";
              success = false;
              break;
            }
          }

          delete[] row_ptrs;
        }

        delete[] rows;
        delete[] cols;
        delete[] vals;

        return success.report(__func__);
      }

    private:
      memory::MemorySpace memspace_;
      MemoryManager       mem_;
    }; // class SparseCooTests
  } // namespace Testing
} // namespace GridKit
