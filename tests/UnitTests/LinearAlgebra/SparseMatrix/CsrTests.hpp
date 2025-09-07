#include <GridKit/LinearAlgebra/SparseMatrix/CSRMatrix.hpp>
#include <GridKit/Utilities/TestHelpers.hpp>
#include <GridKit/Utilities/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    template <class ScalarT, typename IdxT>
    class CsrTests
    {
    public:
      /**
       * @brief Test a simple conversion from COO to CSR
       */
      TestOutcome cooToCsrTest()
      {
        using namespace LinearAlgebra;

        using COO_Matrix = COO_Matrix<ScalarT, IdxT>;
        using CsrMatrix  = CsrMatrix<ScalarT, IdxT>;

        TestStatus success = true;

        std::vector<IdxT>    rows     = {2, 3, 2};
        std::vector<IdxT>    cols     = {3, 3, 2};
        std::vector<ScalarT> vals     = {1, 2, 3};
        const size_t         num_rows = 4;
        const size_t         num_cols = 5;

        COO_Matrix coo(rows, cols, vals, num_rows, num_cols);
        CsrMatrix  csr = CsrMatrix::fromCOO(coo);

        success *= compare(coo, csr);

        return success.report(__func__);
      }

      /**
       * @brief Verifies that you can properly move a COO matrix without copying all of the values
       */
      TestOutcome testCooMove()
      {
        using namespace LinearAlgebra;

        // A dummy struct which simply counts the number of times it has been copied
        struct CopyCounter
        {
          unsigned* counter;

          CopyCounter(unsigned* counter)
            : counter(counter)
          {
          }

          CopyCounter(const CopyCounter& other)
            : counter(other.counter)
          {
            (*counter)++;
          }

          CopyCounter& operator=(const CopyCounter& other)
          {
            counter = other.counter;
            (*counter)++;
            return *this;
          }
        };

        using COO_Matrix = COO_Matrix<CopyCounter, IdxT>;
        using CsrMatrix  = CsrMatrix<CopyCounter, IdxT>;

        TestStatus success = true;

        unsigned counter = 0;

        std::vector<IdxT>        rows = {2, 3, 2};
        std::vector<IdxT>        cols = {3, 3, 2};
        std::vector<CopyCounter> vals(3, CopyCounter(&counter));
        const size_t             num_rows = 4;
        const size_t             num_cols = 5;

        COO_Matrix coo(rows, cols, vals, num_rows, num_cols);
        coo.sortSparse();

        // There should be no copies needed to construct the CSR matrix,
        // since the COO matrix is already sorted.
        counter        = 0;
        CsrMatrix csr  = CsrMatrix::fromCOO(std::move(coo));
        success       *= counter == 0;

        return success.report(__func__);
      }

      /**
       * @brief Test the basic usage of the `LinearAlgebra::CsrBuilder` class with a template matrix
       */
      TestOutcome testCsrBuilderTemplate()
      {
        using namespace LinearAlgebra;

        using COO_Matrix = COO_Matrix<ScalarT, IdxT>;
        using CsrBuilder = CsrBuilder<ScalarT, IdxT>;
        using CsrMatrix  = CsrMatrix<ScalarT, IdxT>;

        TestStatus success = true;

        // Make sure main diagonal elements are added
        std::vector<IdxT>    rows     = {2, 3, 2, 0, 1};
        std::vector<IdxT>    cols     = {3, 3, 2, 0, 1};
        std::vector<ScalarT> vals     = {1, 2, 3, 0, 0};
        const size_t         num_rows = 4;
        const size_t         num_cols = 5;

        COO_Matrix coo(rows, cols, vals, num_rows, num_cols);

        auto builder = CsrBuilder::fromTemplate(CsrMatrix::fromCOO(coo));
        builder.row(2).elem(2, 3).elem(3, 1);
        builder.row(3).elem(3, 2);

        CsrMatrix csr = std::move(builder);

        success *= compare(coo, csr);

        return success.report(__func__);
      }

      /**
       * @brief Test building a CSR matrix the recommended way - from empty at first, then using
       * the original as a template after, with the same building pattern.
       */
      TestOutcome testCsrBuilderComplete()
      {
        using namespace LinearAlgebra;

        using CsrBuilder = CsrBuilder<ScalarT, IdxT>;
        using CsrMatrix  = CsrMatrix<ScalarT, IdxT>;

        TestStatus success = true;

        // Lambda so that we can use both types of builders
        auto build = [](auto builder)
        {
          builder.row(0).elem(2, 3).elem(3, 4).elem(5, 6);
          builder.row(1).elem(0, 1).elem(2, 2);
          builder.row(3).elem(1, 2).elem(3, 3);
          return builder;
        };

        CsrMatrix original_mat = build(CsrBuilder::fromEmpty(4, 6));
        CsrMatrix new_mat      = build(CsrBuilder::fromTemplate(original_mat.clone()));

        success *= compare(original_mat, new_mat);

        return success.report(__func__);
      }

      /**
       * @brief Test building a CSR Matrix with unsorted elements in rows.
       */
      TestOutcome testUnsortedMatrix()
      {
        using namespace LinearAlgebra;

        using CsrBuilder = CsrBuilder<ScalarT, IdxT, true, false>;
        using CsrMatrix  = CsrMatrix<ScalarT, IdxT, false>;

        TestStatus success = true;

        auto build = [](auto builder)
        {
          builder.row(0).elem(2, 3).elem(1, 4).elem(3, 6);
          builder.row(2).elem(3, 1).elem(0, 2).elem(2, 4);

          return builder;
        };

        CsrMatrix original_mat = build(CsrBuilder::fromEmpty(4, 6));
        CsrMatrix new_mat      = build(CsrBuilder::fromTemplate(original_mat.clone()));

        success *= !original_mat.isSorted();
        success *= !new_mat.isSorted();
        success *= compare(original_mat, new_mat);

        return success.report(__func__);
      }

      /**
       * @brief Test the conversion of a matrix from `IS_SORTED == false` to `IS_SORTED == true`.
       */
      TestOutcome testUnsortedToSorted()
      {
        using namespace LinearAlgebra;

        using CsrBuilder = CsrBuilder<ScalarT, IdxT, true, false>;

        TestStatus success = true;

        auto build = [](auto builder)
        {
          builder.row(0).elem(2, 3).elem(1, 4).elem(3, 6);
          builder.row(2).elem(3, 1).elem(0, 2).elem(2, 4);

          return builder;
        };

        CsrMatrix<ScalarT, IdxT, false> original_mat = build(CsrBuilder::fromEmpty(4, 6));

        success *= !original_mat.isSorted();
        success *= !isSorted(original_mat);

        CsrMatrix<ScalarT, IdxT, true> sorted_mat = std::move(original_mat).toSorted();

        success *= sorted_mat.isSorted();
        success *= isSorted(sorted_mat);

        return success.report(__func__);
      }

    private:
      bool compare(LinearAlgebra::CsrMatrix<ScalarT, IdxT>& original_mat, LinearAlgebra::CsrMatrix<ScalarT, IdxT>& new_mat)
      {
        bool comparison = true;

        comparison = comparison
                     && original_mat.numRows() == new_mat.numRows()
                     && original_mat.numCols() == new_mat.numCols()
                     && original_mat.numNonZero() == new_mat.numNonZero();

        for (size_t i = 0; i < original_mat.numRows(); i++)
        {
          for (IdxT j = 0; j < original_mat.numCols(); j++)
          {
            if (original_mat.valueAt(i, j) != new_mat.valueAt(i, j))
            {
              comparison = false;

              std::cerr << "Mismatch at (" << i << "," << j << "): " << original_mat.valueAt(i, j) << " v.s. " << new_mat.valueAt(i, j) << std::endl;
            }
          }
        }

        return comparison;
      }

      bool compare(LinearAlgebra::COO_Matrix<ScalarT, IdxT>& original_mat, LinearAlgebra::CsrMatrix<ScalarT, IdxT>& new_mat)
      {
        bool comparison = true;

        auto [rows, cols, vals] = original_mat.getEntries();

        comparison = comparison
                     && std::get<0>(original_mat.getDimensions()) == new_mat.numRows()
                     && std::get<1>(original_mat.getDimensions()) == new_mat.numCols()
                     && vals.size() == new_mat.numNonZero();

        if (!comparison)
          return comparison;

        ScalarT target;
        // Test all elements - make sure they match
        for (size_t i = 0; i < new_mat.numRows(); i++)
        {
          for (IdxT j = 0; j < new_mat.numCols(); j++)
          {
            target = 0;
            for (size_t k = 0; k < vals.size(); k++)
            {
              if (rows[k] == i && cols[k] == j)
              {
                target = vals[k];
                break;
              }
            }

            comparison = comparison && new_mat.valueAt(i, j) == target;
          }
        }

        return comparison;
      }

      template <bool IS_SORTED>
      bool isSorted(LinearAlgebra::CsrMatrix<ScalarT, IdxT, IS_SORTED>& mat)
      {
        for (size_t i = 1; i < mat.numRows(); i++)
        {
          if (mat.rowIndices()[i + 1] > mat.rowIndices()[i])
          {
            for (IdxT j = mat.rowIndices()[i]; j < mat.rowIndices()[i + 1] - 1; j++)
            {
              // Ascending order of columns
              if (mat.colIndices()[j] >= mat.colIndices()[j + 1])
              {
                return false;
              }
            }
          }
        }

        return true;
      }
    };
  } // namespace Testing
} // namespace GridKit
