#include "LinearAlgebra/SparseMatrix/CSRMatrix.hpp"
#include "Utilities/TestHelpers.hpp"
#include "Utilities/Testing.hpp"

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

        success *= csr.numRows() == num_rows;
        success *= csr.numCols() == num_cols;
        success *= csr.numNonZero() == vals.size();

        ScalarT target;
        // Test all elements - make sure they match
        for (size_t i = 0; i < num_rows; i++)
        {
          for (IdxT j = 0; j < num_cols; j++)
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

            success *= csr.valueAt(i, j) == target;
          }
        }

        return success.report(__func__);
      }
    };
  } // namespace Testing
} // namespace GridKit