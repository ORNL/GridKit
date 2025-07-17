#pragma once

namespace GridKit
{
  namespace LinearAlgebra
  {
    template <typename ScalarT, typename IdxT>
    class CsrMatrixRef
    {
    protected:
      IdxT n_rows_;    ///< The number of rows in this matrix
      IdxT n_columns_; ///< The number of columns in this matrix
      IdxT nnz_;       ///< The number of nonzeroes in this matrix

      IdxT*    row_pointers_;   ///< The row pointers of the CSR matrix
      IdxT*    column_indices_; ///< The column indices of the CSR matrix
      ScalarT* values_;         ///< The values of the CSR matrix
    };

    template <typename ScalarT, typename IdxT>
    class CsrMatrix : public CsrMatrixRef<ScalarT, IdxT>
    {
    protected:
    public:
      // ok
    };
  } // namespace LinearAlgebra
} // namespace GridKit
