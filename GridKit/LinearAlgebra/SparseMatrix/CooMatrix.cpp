
#include "CooMatrix.hpp"

#include <algorithm>
#include <cassert>
#include <iomanip>
#include <limits>

namespace GridKit
{
  namespace LinearAlgebra
  {
    template <typename RealT, typename IdxT>
    CooMatrix<RealT, IdxT>::CooMatrix()
    {
    }

    /**
     * @brief basic constructor. It DOES NOT allocate any memory!
     *
     * @param[in] n   - number of rows
     * @param[in] m   - number of columns
     * @param[in] nnz_ - number of non-zeros
     */
    template <typename RealT, typename IdxT>
    CooMatrix<RealT, IdxT>::CooMatrix(IdxT n,
                                      IdxT m,
                                      IdxT nnz)
      : n_{n},
        m_{m},
        nnz_{nnz}
    {
      setNotUpdated();

      // set everything to nullptr
      h_row_data_ = nullptr;
      h_col_data_ = nullptr;
      h_val_data_ = nullptr;

      d_row_data_ = nullptr;
      d_col_data_ = nullptr;
      d_val_data_ = nullptr;

      owns_cpu_sparsity_pattern_ = false;
      owns_cpu_values_           = false;

      owns_gpu_sparsity_pattern_ = false;
      owns_gpu_values_           = false;
    }

    /**
     * @brief Hijacking constructor
     *
     * @param[in] n
     * @param[in] m
     * @param[in] nnz_
     * @param[in,out] rows
     * @param[in,out] cols
     * @param[in,out] vals
     * @param[in] memspace_src
     * @param[in] memspace_dst
     */
    template <typename RealT, typename IdxT>
    CooMatrix<RealT, IdxT>::CooMatrix(IdxT                n,
                                      IdxT                m,
                                      IdxT                nnz_,
                                      IdxT**              rows,
                                      IdxT**              cols,
                                      RealT**             vals,
                                      memory::MemorySpace memspace_src,
                                      memory::MemorySpace memspace_dst)
      : CooMatrix(n, m, nnz_)
    {
      int control = -1;
      if ((memspace_src == memory::HOST) && (memspace_dst == memory::HOST))
      {
        control = 0;
      }
      if ((memspace_src == memory::HOST) && (memspace_dst == memory::DEVICE))
      {
        control = 1;
      }
      if ((memspace_src == memory::DEVICE) && (memspace_dst == memory::HOST))
      {
        control = 2;
      }
      if ((memspace_src == memory::DEVICE) && (memspace_dst == memory::DEVICE))
      {
        control = 3;
      }

      switch (control)
      {
      case 0: // cpu->cpu
        // Set host data
        h_row_data_                = *rows;
        h_col_data_                = *cols;
        h_val_data_                = *vals;
        h_data_updated_            = true;
        owns_cpu_sparsity_pattern_ = true;
        owns_cpu_values_           = true;
        // Set device data to null
        if (d_row_data_ || d_col_data_ || d_val_data_)
        {
          std::cerr << "Device data unexpectedly allocated. "
                    << "Possible bug in matrix::Sparse class.\n";
        }
        d_row_data_                = nullptr;
        d_col_data_                = nullptr;
        d_val_data_                = nullptr;
        d_data_updated_            = false;
        owns_gpu_sparsity_pattern_ = false;
        owns_gpu_values_           = false;
        // Hijack data from the source
        *rows                      = nullptr;
        *cols                      = nullptr;
        *vals                      = nullptr;
        break;
      case 2: // gpu->cpu
        // Set device data and copy it to host
        d_row_data_                = *rows;
        d_col_data_                = *cols;
        d_val_data_                = *vals;
        d_data_updated_            = true;
        owns_gpu_sparsity_pattern_ = true;
        owns_gpu_values_           = true;
        syncData(memspace_dst);
        // Hijack data from the source
        *rows = nullptr;
        *cols = nullptr;
        *vals = nullptr;
        break;
      case 1: // cpu->gpu
        // Set host data and copy it to device
        h_row_data_                = *rows;
        h_col_data_                = *cols;
        h_val_data_                = *vals;
        h_data_updated_            = true;
        owns_cpu_sparsity_pattern_ = true;
        owns_cpu_values_           = true;
        syncData(memspace_dst);

        // Hijack data from the source
        *rows = nullptr;
        *cols = nullptr;
        *vals = nullptr;
        break;
      case 3: // gpu->gpu
        // Set device data
        d_row_data_                = *rows;
        d_col_data_                = *cols;
        d_val_data_                = *vals;
        d_data_updated_            = true;
        owns_gpu_sparsity_pattern_ = true;
        owns_gpu_values_           = true;
        // Set host data to null
        if (h_row_data_ || h_col_data_ || h_val_data_)
        {
          std::cerr << "Host data unexpectedly allocated. "
                    << "Possible bug in matrix::Sparse class.\n";
        }
        h_row_data_                = nullptr;
        h_col_data_                = nullptr;
        h_val_data_                = nullptr;
        h_data_updated_            = false;
        owns_cpu_sparsity_pattern_ = false;
        owns_cpu_values_           = false;
        // Hijack data from the source
        *rows                      = nullptr;
        *cols                      = nullptr;
        *vals                      = nullptr;
        break;
      default:
        std::cerr << "CooMatrix constructor failed! "
                  << "Possible bug in memory spaces setting.\n";
        break;
      }
    }

    template <typename RealT, typename IdxT>
    CooMatrix<RealT, IdxT>::~CooMatrix()
    {
      this->destroyMatrixData(memory::HOST);
      this->destroyMatrixData(memory::DEVICE);
      delete[] map_to_sorted_;
      delete[] map_to_dedup_;
    }

    /**
     * @brief set the matrix update flags to false (for both HOST and DEVICE).
     */
    template <typename RealT, typename IdxT>
    void CooMatrix<RealT, IdxT>::setNotUpdated()
    {
      h_data_updated_ = false;
      d_data_updated_ = false;
    }

    /**
     * @brief get number of matrix rows
     *
     * @return number of matrix rows.
     */
    template <typename RealT, typename IdxT>
    IdxT CooMatrix<RealT, IdxT>::getNumRows() const
    {
      return this->n_;
    }

    /**
     * @brief get number of matrix columns
     *
     * @return number of matrix columns.
     */
    template <typename RealT, typename IdxT>
    IdxT CooMatrix<RealT, IdxT>::getNumColumns() const
    {
      return this->m_;
    }

    /**
     * @brief get number of non-zeros in the matrix.
     *
     * @return number of non-zeros.
     */
    template <typename RealT, typename IdxT>
    IdxT CooMatrix<RealT, IdxT>::getNnz() const
    {
      return this->nnz_;
    }

    /**
     * @brief Set the pointers for matrix row, column, value data.
     *
     * Useful if interfacing with other codes - this function only assigns
     * pointers, but it does not allocate nor copy anything. The data ownership
     * flags would be set to false (default).
     *
     * @param[in] row_data - pointer to row data (array of integers)
     * @param[in] col_data - pointer to column data (array of integers)
     * @param[in] val_data - pointer to value data (array of real numbers)
     * @param[in] memspace - memory space (HOST or DEVICE) of incoming data
     *
     * @return 0 if successful, 1 if not.
     */
    template <typename RealT, typename IdxT>
    int CooMatrix<RealT, IdxT>::setDataPointers(IdxT*               row_data,
                                                IdxT*               col_data,
                                                RealT*              val_data,
                                                memory::MemorySpace memspace)
    {
      using namespace memory;

      setNotUpdated();

      switch (memspace)
      {
      case HOST:
        if (owns_cpu_sparsity_pattern_ && (h_row_data_ || h_col_data_))
        {
          std::cerr << "Trying to set matrix host data, but the data already set!\n";
          std::cerr << "Ignoring setDataPointers function call ...\n";
          return 1;
        }
        if (owns_cpu_values_ && h_val_data_)
        {
          std::cerr << "Trying to set matrix host values, but the values already set!\n";
          std::cerr << "Ignoring setValuesPointer function call ...\n";
          return 1;
        }
        h_row_data_                = row_data;
        h_col_data_                = col_data;
        h_val_data_                = val_data;
        h_data_updated_            = true;
        owns_cpu_sparsity_pattern_ = false;
        owns_cpu_values_           = false;
        break;
      case DEVICE:
        if (owns_gpu_sparsity_pattern_ && (d_row_data_ || d_col_data_))
        {
          std::cerr << "Trying to set matrix host data, but the data already set!\n";
          std::cerr << "Ignoring setDataPointers function call ...\n";
          return 1;
        }
        if (owns_gpu_values_ && d_val_data_)
        {
          std::cerr << "Trying to set matrix device values, but the values already set!\n";
          std::cerr << "Ignoring setValuesPointer function call ...\n";
          return 1;
        }
        d_row_data_                = row_data;
        d_col_data_                = col_data;
        d_val_data_                = val_data;
        d_data_updated_            = true;
        owns_gpu_sparsity_pattern_ = false;
        owns_gpu_values_           = false;
        break;
      }
      return 0;
    }

    /**
     * @brief destroy matrix data (HOST or DEVICE) if the matrix owns it
     * (will attempt to destroy all three arrays).
     *
     * @param[in] memspace - memory space (HOST or DEVICE) of incoming data
     *
     * @return 0 if successful, -1 if not.
     *
     */
    template <typename RealT, typename IdxT>
    int CooMatrix<RealT, IdxT>::destroyMatrixData(memory::MemorySpace memspace)
    {
      using namespace memory;
      switch (memspace)
      {
      case HOST:
        if (owns_cpu_sparsity_pattern_)
        {
          delete[] h_row_data_;
          delete[] h_col_data_;
          h_row_data_ = nullptr;
          h_col_data_ = nullptr;
        }
        if (owns_cpu_values_)
        {
          delete[] h_val_data_;
          h_val_data_ = nullptr;
        }
        return 0;
      case DEVICE:
        if (owns_gpu_sparsity_pattern_)
        {
          mem_.deleteOnDevice(d_row_data_);
          mem_.deleteOnDevice(d_col_data_);
          d_row_data_ = nullptr;
          d_col_data_ = nullptr;
        }
        if (owns_gpu_values_)
        {
          mem_.deleteOnDevice(d_val_data_);
          d_val_data_ = nullptr;
        }
        return 0;
      default:
        return -1;
      }
    }

    /**
     * @brief Sort COO entries by (row, col), merge duplicate
     *        entries, and build CSR row pointers.
     *
     * @note This function operates on host (CPU) data only.
     *
     * @return Pointer to CSR row pointer array
     */
    template <typename RealT, typename IdxT>
    IdxT* CooMatrix<RealT, IdxT>::getCsrRowData()
    {
      if (!h_data_updated_)
      {
        std::cerr << "CooMatrix::getCsrRowData requires up-to-date host data, but host is out of date!\n";
        assert(h_data_updated_);
        return nullptr;
      }

      map_to_sorted_ = new IdxT[nnz_];

      for (IdxT i = 0; i < nnz_; i++)
      {
        map_to_sorted_[i] = i;
      }

      // Sort by (row, col)
      std::sort(map_to_sorted_, map_to_sorted_ + nnz_, [this](IdxT a, IdxT b)
                { if (h_row_data_[a] != h_row_data_[b]) return h_row_data_[a] < h_row_data_[b];
                  return h_col_data_[a] < h_col_data_[b]; });

      IdxT*  row_data_tmp = new IdxT[nnz_];
      IdxT*  col_data_tmp = new IdxT[nnz_];
      RealT* val_data_tmp = new RealT[nnz_];

      for (IdxT i = 0; i < nnz_; i++)
      {
        row_data_tmp[i] = h_row_data_[map_to_sorted_[i]];
        col_data_tmp[i] = h_col_data_[map_to_sorted_[i]];
        val_data_tmp[i] = h_val_data_[map_to_sorted_[i]];
      }
      for (IdxT i = 0; i < nnz_; i++)
      {
        h_row_data_[i] = row_data_tmp[i];
        h_col_data_[i] = col_data_tmp[i];
        h_val_data_[i] = val_data_tmp[i];
      }
      delete[] row_data_tmp;
      delete[] col_data_tmp;
      delete[] val_data_tmp;

      map_to_dedup_      = new IdxT[nnz_];
      IdxT* csr_row_data = new IdxT[n_ + 1]();

      map_to_dedup_[0] = 0;
      csr_row_data[h_row_data_[0] + 1]++;

      // Deduplicate sroted entries
      IdxT w = 0;
      for (IdxT i = 1; i < nnz_; i++)
      {
        if (h_row_data_[i] == h_row_data_[w] && h_col_data_[i] == h_col_data_[w])
        {
          h_val_data_[w]   += h_val_data_[i];
          map_to_dedup_[i]  = w;
        }
        else
        {
          ++w;
          h_row_data_[w]   = h_row_data_[i];
          h_col_data_[w]   = h_col_data_[i];
          h_val_data_[w]   = h_val_data_[i];
          map_to_dedup_[i] = w;
          csr_row_data[h_row_data_[w] + 1]++;
        }
      }
      // nnz after deduplication
      nnz_ = w + 1;

      for (IdxT i = 0; i < n_; i++)
      {
        csr_row_data[i + 1] += csr_row_data[i];
      }

      return csr_row_data;
    }

    template <typename RealT, typename IdxT>
    const IdxT* CooMatrix<RealT, IdxT>::getMapToSorted() const
    {
      return map_to_sorted_;
    }

    template <typename RealT, typename IdxT>
    const IdxT* CooMatrix<RealT, IdxT>::getMapToDeduplicated() const
    {
      return map_to_dedup_;
    }

    template <typename RealT, typename IdxT>
    IdxT* CooMatrix<RealT, IdxT>::getRowData(memory::MemorySpace memspace)
    {
      using namespace memory;

      switch (memspace)
      {
      case HOST:
        return this->h_row_data_;
      case DEVICE:
        return this->d_row_data_;
      default:
        return nullptr;
      }
    }

    template <typename RealT, typename IdxT>
    IdxT* CooMatrix<RealT, IdxT>::getColData(memory::MemorySpace memspace)
    {
      using namespace memory;

      switch (memspace)
      {
      case HOST:
        return this->h_col_data_;
      case DEVICE:
        return this->d_col_data_;
      default:
        return nullptr;
      }
    }

    template <typename RealT, typename IdxT>
    RealT* CooMatrix<RealT, IdxT>::getValues(memory::MemorySpace memspace)
    {
      using namespace memory;

      switch (memspace)
      {
      case HOST:
        return this->h_val_data_;
      case DEVICE:
        return this->d_val_data_;
      default:
        return nullptr;
      }
    }

    /**
     * @brief Sync data in memspace with the updated memory space.
     *
     * @param memspace - memory space to be synced up (HOST or DEVICE)
     * @return int - 0 if successful, error code otherwise
     *
     * @pre The memory space other than `memspace` must be up-to-date. Otherwise,
     * this function will return an error.
     *
     * @see CooMatrix::setUpdated
     */
    template <typename RealT, typename IdxT>
    int CooMatrix<RealT, IdxT>::syncData(memory::MemorySpace memspace)
    {
      using namespace memory;

      switch (memspace)
      {
      case HOST:
        // check if we need to copy or not
        assert(((h_row_data_ == nullptr) == (h_col_data_ == nullptr)) && "In CooMatrix::syncData one of host row or column data is null!\n");

        if (h_data_updated_)
        {
          std::cerr << "CooMatrix::syncData is trying to sync host, but host already up to date!\n";
          assert(!h_data_updated_);
          return 1;
        }
        if (!d_data_updated_)
        {
          std::cerr << "CooMatrix::syncData is trying to sync host with device, but device is out of date!\n"
                    << "See CooMatrix::syncData documentation\n.";
          assert(d_data_updated_);
        }
        if ((h_row_data_ == nullptr) && (h_col_data_ == nullptr))
        {
          h_row_data_                = new IdxT[static_cast<size_t>(n_ + 1)];
          h_col_data_                = new IdxT[static_cast<size_t>(nnz_)];
          owns_cpu_sparsity_pattern_ = true;
        }
        if (h_val_data_ == nullptr)
        {
          h_val_data_      = new RealT[static_cast<size_t>(nnz_)];
          owns_cpu_values_ = true;
        }
        mem_.copyArrayDeviceToHost(h_row_data_, d_row_data_, n_ + 1);
        mem_.copyArrayDeviceToHost(h_col_data_, d_col_data_, nnz_);
        mem_.copyArrayDeviceToHost(h_val_data_, d_val_data_, nnz_);
        h_data_updated_ = true;
        return 0;
      case DEVICE:
        assert(((d_row_data_ == nullptr) == (d_col_data_ == nullptr)) && "In CooMatrix::syncData one of device row or column data is null!\n");

        if (d_data_updated_)
        {
          std::cerr << "CooMatrix::syncData is trying to sync device, but device already up to date!\n";
          assert(!d_data_updated_);
          return 1;
        }
        if (!h_data_updated_)
        {
          std::cerr << "CooMatrix::syncData is trying to sync device with host, but host is out of date!\n"
                    << "See CooMatrix::syncData documentation\n.";
          assert(h_data_updated_);
        }
        if ((d_row_data_ == nullptr) && (d_col_data_ == nullptr))
        {
          mem_.allocateArrayOnDevice(&d_row_data_, n_ + 1);
          mem_.allocateArrayOnDevice(&d_col_data_, nnz_);
          owns_gpu_sparsity_pattern_ = true;
        }
        if (d_val_data_ == nullptr)
        {
          mem_.allocateArrayOnDevice(&d_val_data_, nnz_);
          owns_gpu_values_ = true;
        }
        mem_.copyArrayHostToDevice(d_row_data_, h_row_data_, n_ + 1);
        mem_.copyArrayHostToDevice(d_col_data_, h_col_data_, nnz_);
        mem_.copyArrayHostToDevice(d_val_data_, h_val_data_, nnz_);
        d_data_updated_ = true;
        return 0;
      default:
        return 1;
      } // switch
    }

    /**
     * @brief Prints matrix data.
     *
     * @param out - Output stream where the matrix data is printed
     */
    template <typename RealT, typename IdxT>
    void CooMatrix<RealT, IdxT>::print(std::ostream& out, IdxT indexing_base)
    {
      out << std::scientific << std::setprecision(std::numeric_limits<RealT>::digits10);
      for (IdxT i = 0; i < n_; ++i)
      {
        for (IdxT j = h_row_data_[i]; j < h_row_data_[i + 1]; ++j)
        {
          out << i + indexing_base << " "
              << h_col_data_[j] + indexing_base << " "
              << h_val_data_[j] << "\n";
        }
      }
    }

    template class CooMatrix<double, long int>;
    template class CooMatrix<double, size_t>;
    template class CooMatrix<double, int>;
  } // namespace LinearAlgebra
} // namespace GridKit