#pragma once

#include <cstdint>
#include <cstring>
#include <iostream>

#include <GridKit/LinearAlgebra/MemoryUtils.hpp>

namespace GridKit
{
  namespace LinearAlgebra
  {

    template <typename RealT, typename IdxT>
    class CsrMatrix
    {
    public:
      CsrMatrix();

      CsrMatrix(IdxT n, IdxT m, IdxT nnz);

      CsrMatrix(IdxT                n,
                IdxT                m,
                IdxT                nnz,
                IdxT**              rows,
                IdxT**              cols,
                RealT**             vals,
                memory::MemorySpace memspace_src = memory::HOST,
                memory::MemorySpace memspace_dst = memory::HOST);

      ~CsrMatrix();

      // accessors
      IdxT getNumRows();
      IdxT getNumColumns();
      IdxT getNnz();

      void setNnz(IdxT nnz_new); // for resetting when removing duplicates
      int  setUpdated(memory::MemorySpace what);

      IdxT*  getRowData(memory::MemorySpace memspace = memory::HOST);
      IdxT*  getColData(memory::MemorySpace memspace = memory::HOST);
      RealT* getValues(memory::MemorySpace memspace = memory::HOST);

      int copyDataFrom(const IdxT*         row_data,
                       const IdxT*         col_data,
                       const RealT*        val_data,
                       memory::MemorySpace memspace_in,
                       memory::MemorySpace memspace_out);
      int copyDataFrom(const IdxT*         row_data,
                       const IdxT*         col_data,
                       const RealT*        val_data,
                       IdxT                new_nnz,
                       memory::MemorySpace memspace_in,
                       memory::MemorySpace memspace_out);

      int allocateMatrixData(memory::MemorySpace memspace);
      int setDataPointers(IdxT*               row_data,
                          IdxT*               col_data,
                          RealT*              val_data,
                          memory::MemorySpace memspace);

      int destroyMatrixData(memory::MemorySpace memspace);

      void print(std::ostream& file_out = std::cout, IdxT indexing_base = 0);

      int syncData(memory::MemorySpace memspace_out);

      // update Values just updates values; it allocates if necessary.
      // values have the same dimensions between different formats
      virtual int copyValues(const RealT*        new_vals,
                             memory::MemorySpace memspace_in,
                             memory::MemorySpace memspace_out);

      // set new values just sets the pointer, use caution.
      virtual int setValuesPointer(RealT*              new_vals,
                                   memory::MemorySpace memspace);

    private:
      IdxT n_{0};   ///< number of rows
      IdxT m_{0};   ///< number of columns
      IdxT nnz_{0}; ///< number of non-zeros

      // host data
      IdxT*  h_row_data_{nullptr};   ///< row data (HOST)
      IdxT*  h_col_data_{nullptr};   ///< column data (HOST)
      RealT* h_val_data_{nullptr};   ///< value data (HOST)
      bool   h_data_updated_{false}; ///< HOST update flag

      // gpu data
      IdxT*  d_row_data_{nullptr};   ///< row data (DEVICE)
      IdxT*  d_col_data_{nullptr};   ///< column data (DEVICE)
      RealT* d_val_data_{nullptr};   ///< value data (DEVICE)
      bool   d_data_updated_{false}; ///< DEVICE update flag

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
