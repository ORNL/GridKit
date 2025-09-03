#include <GridKit/LinearAlgebra/SparseMatrix/CSRMatrix.hpp>
#include <GridKit/Utilities/TestHelpers.hpp>
#include <GridKit/Utilities/Testing.hpp>

namespace GridKit
{
  namespace Testing
  {
    template <class ScalarT, typename IdxT>
    class CSRTests
    {
    public:
      /**
       * @brief Test a simple conversion from COO to CSR
       */
      TestOutcome cooToCsrTest()
      {
        using namespace LinearAlgebra;

        using COO_Matrix = COO_Matrix<ScalarT, IdxT>;
        using CSRMatrix  = CSRMatrix<ScalarT, IdxT>;

        TestStatus success = true;

        std::vector<IdxT>    rows     = {2, 3, 2};
        std::vector<IdxT>    cols     = {3, 3, 2};
        std::vector<ScalarT> vals     = {1, 2, 3};
        const size_t         num_rows = 4;
        const size_t         num_cols = 5;

        COO_Matrix coo(rows, cols, vals, num_rows, num_cols);
        CSRMatrix  csr = CSRMatrix::fromCOO(coo);

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

          CopyCounter(unsigned* counter) : counter(counter)
          {
          }

          CopyCounter(const CopyCounter& other) : counter(other.counter)
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
        using CSRMatrix  = CSRMatrix<CopyCounter, IdxT>;

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
        CSRMatrix csr  = CSRMatrix::fromCOO(std::move(coo));
        success       *= counter == 0;

        return success.report(__func__);
      }

      /**
       * @brief Test the basic usage of the `LinearAlgebra::CSRBuilder` class with a template matrix
       */
      TestOutcome testCsrBuilderTemplate()
      {
        using namespace LinearAlgebra;

        using COO_Matrix = COO_Matrix<ScalarT, IdxT>;
        using CSRBuilder = CSRBuilder<ScalarT, IdxT>;
        using CSRMatrix  = CSRMatrix<ScalarT, IdxT>;

        TestStatus success = true;

        // Make sure main diagonal elements are added
        std::vector<IdxT>    rows     = {2, 3, 2, 0, 1};
        std::vector<IdxT>    cols     = {3, 3, 2, 0, 1};
        std::vector<ScalarT> vals     = {1, 2, 3, 0, 0};
        const size_t         num_rows = 4;
        const size_t         num_cols = 5;

        COO_Matrix coo(rows, cols, vals, num_rows, num_cols);

        auto builder = CSRBuilder::fromTemplate(CSRMatrix::fromCOO(coo));
        builder.row(2).elem(2, 3).elem(3, 1);
        builder.row(3).elem(3, 2);

        CSRMatrix csr = std::move(builder);

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

        using CSRBuilder = CSRBuilder<ScalarT, IdxT>;
        using CSRMatrix  = CSRMatrix<ScalarT, IdxT>;

        TestStatus success = true;

        // Lambda so that we can use both types of builders
        auto build = [](auto builder)
        {
          builder.row(0).elem(2, 3).elem(3, 4).elem(5, 6);
          builder.row(1).elem(0, 1).elem(2, 2);
          builder.row(3).elem(1, 2).elem(3, 3);
          return builder;
        };

        CSRMatrix original_mat = build(CSRBuilder::fromEmpty(4, 6));
        CSRMatrix new_mat      = build(CSRBuilder::fromTemplate(original_mat.clone()));

        success *= compare(original_mat, new_mat);

        return success.report(__func__);
      }

    private:
      bool compare(LinearAlgebra::CSRMatrix<ScalarT, IdxT>& original_mat, LinearAlgebra::CSRMatrix<ScalarT, IdxT>& new_mat)
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
            comparison = comparison && original_mat.valueAt(i, j) == new_mat.valueAt(i, j);
          }
        }

        return comparison;
      }

      bool compare(LinearAlgebra::COO_Matrix<ScalarT, IdxT>& original_mat, LinearAlgebra::CSRMatrix<ScalarT, IdxT>& new_mat)
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
    };
  } // namespace Testing
} // namespace GridKit