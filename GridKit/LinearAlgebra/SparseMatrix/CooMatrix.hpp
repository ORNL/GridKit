#pragma once

#include <cstdint>
#include <cstring>
#include <iostream>

#include <GridKit/MemoryUtilities/MemoryUtils.hpp>

namespace GridKit
{
  namespace LinearAlgebra
  {

    template <typename RealT, typename IdxT>
    class CooMatrix
    {
    public:
      CooMatrix();

      CooMatrix(IdxT n, IdxT m, IdxT nnz);

      CooMatrix(IdxT                n,
                IdxT                m,
                IdxT                nnz,
                IdxT**              rows,
                IdxT**              cols,
                RealT**             vals,
                memory::MemorySpace memspace_src = memory::HOST,
                memory::MemorySpace memspace_dst = memory::HOST);

      ~CooMatrix();

      // accessors
      IdxT getNumRows() const;
      IdxT getNumColumns() const;
      IdxT getNnz() const;

      const IdxT* getMapToSorted() const;
      const IdxT* getMapToDeduplicated() const;

      IdxT*  getRowData(memory::MemorySpace memspace = memory::HOST);
      IdxT*  getColData(memory::MemorySpace memspace = memory::HOST);
      RealT* getValues(memory::MemorySpace memspace = memory::HOST);

      int setDataPointers(IdxT*               row_data,
                          IdxT*               col_data,
                          RealT*              val_data,
                          memory::MemorySpace memspace);

      int destroyMatrixData(memory::MemorySpace memspace);

      void print(std::ostream& file_out = std::cout, IdxT indexing_base = 0);

      int syncData(memory::MemorySpace memspace_out);

      IdxT* getCsrRowData();

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

      IdxT* map_to_sorted_ = {nullptr}; ///< map from orginal to sorted
      IdxT* map_to_dedup_  = {nullptr}; ///< map from sorted to deduplicated

      MemoryManager mem_; ///< Device memory manager object
    };
  } // namespace LinearAlgebra
} // namespace GridKit
