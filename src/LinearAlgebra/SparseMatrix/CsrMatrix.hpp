#pragma once

#include <cstdint>
#include <cstring>
#include <iostream>

#include <LinearAlgebra/MemoryUtils.hpp>

namespace GridKit
{
  namespace LinearAlgebra
  {
    using real_type  = double;
    using index_type = std::int32_t;

    class CsrMatrix
    {
    public:
      CsrMatrix();

      CsrMatrix(index_type n, index_type m, index_type nnz);

      CsrMatrix(index_type          n,
                index_type          m,
                index_type          nnz,
                index_type**        rows,
                index_type**        cols,
                real_type**         vals,
                memory::MemorySpace memspaceSrc = memory::HOST,
                memory::MemorySpace memspaceDst = memory::HOST);

      ~CsrMatrix();

      // accessors
      index_type getNumRows();
      index_type getNumColumns();
      index_type getNnz();

      void setNnz(index_type nnz_new); // for resetting when removing duplicates
      int  setUpdated(memory::MemorySpace what);

      index_type* getRowData(memory::MemorySpace memspace = memory::HOST);
      index_type* getColData(memory::MemorySpace memspace = memory::HOST);
      real_type*  getValues(memory::MemorySpace memspace = memory::HOST);

      int copyDataFrom(const index_type*   row_data,
                       const index_type*   col_data,
                       const real_type*    val_data,
                       memory::MemorySpace memspaceIn,
                       memory::MemorySpace memspaceOut);
      int copyDataFrom(const index_type*   row_data,
                       const index_type*   col_data,
                       const real_type*    val_data,
                       index_type          new_nnz,
                       memory::MemorySpace memspaceIn,
                       memory::MemorySpace memspaceOut);

      int allocateMatrixData(memory::MemorySpace memspace);
      int setDataPointers(index_type*         row_data,
                          index_type*         col_data,
                          real_type*          val_data,
                          memory::MemorySpace memspace);

      int destroyMatrixData(memory::MemorySpace memspace);

      void print(std::ostream& file_out = std::cout, index_type indexing_base = 0);

      int syncData(memory::MemorySpace memspaceOut);

      // update Values just updates values; it allocates if necessary.
      // values have the same dimensions between different formats
      virtual int copyValues(const real_type*    new_vals,
                             memory::MemorySpace memspaceIn,
                             memory::MemorySpace memspaceOut);

      // set new values just sets the pointer, use caution.
      virtual int setValuesPointer(real_type*          new_vals,
                                   memory::MemorySpace memspace);

    private:
      index_type n_{0};   ///< number of rows
      index_type m_{0};   ///< number of columns
      index_type nnz_{0}; ///< number of non-zeros

      // host data
      index_type* h_row_data_{nullptr};   ///< row data (HOST)
      index_type* h_col_data_{nullptr};   ///< column data (HOST)
      real_type*  h_val_data_{nullptr};   ///< value data (HOST)
      bool        h_data_updated_{false}; ///< HOST update flag

      // gpu data
      index_type* d_row_data_{nullptr};   ///< row data (DEVICE)
      index_type* d_col_data_{nullptr};   ///< column data (DEVICE)
      real_type*  d_val_data_{nullptr};   ///< value data (DEVICE)
      bool        d_data_updated_{false}; ///< DEVICE update flag

      void setNotUpdated();

      // Data ownership flags
      bool owns_cpu_sparsity_pattern_{false}; ///< for row/col data
      bool owns_cpu_values_{false};           ///< for nonzero values

      bool owns_gpu_sparsity_pattern_{false}; ///< for row/col data
      bool owns_gpu_values_{false};           ///< for nonzero values

      MemoryHandler mem_; ///< Device memory manager object
    };
  } // namespace LinearAlgebra
} // namespace GridKit
