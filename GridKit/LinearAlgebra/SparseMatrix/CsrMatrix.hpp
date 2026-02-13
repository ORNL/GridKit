#pragma once

#include <cassert>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>

#include <GridKit/LinearAlgebra/MemoryUtils.hpp>

namespace GridKit
{
  namespace LinearAlgebra
  {

    template <typename RealT, typename IdxT>
    class CsrMatrix
    {
    public:
      CsrMatrix()
      {
      }

      /**
       * @brief basic constructor. It DOES NOT allocate any memory!
       *
       * @param[in] n   - number of rows
       * @param[in] m   - number of columns
       * @param[in] nnz - number of non-zeros
       */
      CsrMatrix(IdxT n, IdxT m, IdxT nnz)
        : n_{n},
          m_{m},
          nnz_{nnz}
      {
        setNotUpdated();

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
       * @param[in] nnz
       * @param[in,out] rows
       * @param[in,out] cols
       * @param[in,out] vals
       * @param[in] memspace_src
       * @param[in] memspace_dst
       */
      CsrMatrix(IdxT                n,
                IdxT                m,
                IdxT                nnz,
                IdxT**              rows,
                IdxT**              cols,
                RealT**             vals,
                memory::MemorySpace memspace_src = memory::HOST,
                memory::MemorySpace memspace_dst = memory::HOST)
        : CsrMatrix(n, m, nnz)
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
          h_row_data_                = *rows;
          h_col_data_                = *cols;
          h_val_data_                = *vals;
          h_data_updated_            = true;
          owns_cpu_sparsity_pattern_ = true;
          owns_cpu_values_           = true;
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
          *rows                      = nullptr;
          *cols                      = nullptr;
          *vals                      = nullptr;
          break;
        case 2: // gpu->cpu
          d_row_data_                = *rows;
          d_col_data_                = *cols;
          d_val_data_                = *vals;
          d_data_updated_            = true;
          owns_gpu_sparsity_pattern_ = true;
          owns_gpu_values_           = true;
          syncData(memspace_dst);
          *rows = nullptr;
          *cols = nullptr;
          *vals = nullptr;
          break;
        case 1: // cpu->gpu
          h_row_data_                = *rows;
          h_col_data_                = *cols;
          h_val_data_                = *vals;
          h_data_updated_            = true;
          owns_cpu_sparsity_pattern_ = true;
          owns_cpu_values_           = true;
          syncData(memspace_dst);
          *rows = nullptr;
          *cols = nullptr;
          *vals = nullptr;
          break;
        case 3: // gpu->gpu
          d_row_data_                = *rows;
          d_col_data_                = *cols;
          d_val_data_                = *vals;
          d_data_updated_            = true;
          owns_gpu_sparsity_pattern_ = true;
          owns_gpu_values_           = true;
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
          *rows                      = nullptr;
          *cols                      = nullptr;
          *vals                      = nullptr;
          break;
        default:
          std::cerr << "CsrMatrix constructor failed! "
                    << "Possible bug in memory spaces setting.\n";
          break;
        }
      }

      ~CsrMatrix()
      {
        this->destroyMatrixData(memory::HOST);
        this->destroyMatrixData(memory::DEVICE);
      }

      // accessors

      IdxT getNumRows() const
      {
        return n_;
      }

      IdxT getNumColumns() const
      {
        return m_;
      }

      IdxT getNnz() const
      {
        return nnz_;
      }

      IdxT getRows() const
      {
        return n_;
      }

      IdxT getCols() const
      {
        return m_;
      }

      const IdxT* getRowPtr() const
      {
        return h_row_data_;
      }

      const IdxT* getColInd() const
      {
        return h_col_data_;
      }

      const RealT* getValues() const
      {
        return h_val_data_;
      }

      void setNnz(IdxT nnz_new)
      {
        nnz_ = nnz_new;
      }

      /**
       * @brief Tags `memspace` as updated.
       *
       * @param[in] memspace - memory space (HOST or DEVICE) to set to "updated"
       *
       * @return 0 if successful, -1 if not.
       *
       * @warning This is an expert-level function. Use only if you know what you
       * are doing.
       */
      int setUpdated(memory::MemorySpace memspace)
      {
        using namespace memory;
        switch (memspace)
        {
        case HOST:
          h_data_updated_ = true;
          d_data_updated_ = false;
          break;
        case DEVICE:
          d_data_updated_ = true;
          h_data_updated_ = false;
          break;
        }
        return 0;
      }

      IdxT* getRowData(memory::MemorySpace memspace = memory::HOST)
      {
        using namespace memory;
        switch (memspace)
        {
        case HOST:
          return h_row_data_;
        case DEVICE:
          return d_row_data_;
        default:
          return nullptr;
        }
      }

      IdxT* getColData(memory::MemorySpace memspace = memory::HOST)
      {
        using namespace memory;
        switch (memspace)
        {
        case HOST:
          return h_col_data_;
        case DEVICE:
          return d_col_data_;
        default:
          return nullptr;
        }
      }

      RealT* getValues(memory::MemorySpace memspace = memory::HOST)
      {
        using namespace memory;
        switch (memspace)
        {
        case HOST:
          return h_val_data_;
        case DEVICE:
          return d_val_data_;
        default:
          return nullptr;
        }
      }

      int copyDataFrom(const IdxT*         row_data,
                       const IdxT*         col_data,
                       const RealT*        val_data,
                       memory::MemorySpace memspace_in,
                       memory::MemorySpace memspace_out)
      {
        IdxT nnz_current = nnz_;
        setNotUpdated();
        int control = -1;
        if ((memspace_in == memory::HOST) && (memspace_out == memory::HOST))
        {
          control = 0;
        }
        if ((memspace_in == memory::HOST) && ((memspace_out == memory::DEVICE)))
        {
          control = 1;
        }
        if (((memspace_in == memory::DEVICE)) && (memspace_out == memory::HOST))
        {
          control = 2;
        }
        if (((memspace_in == memory::DEVICE)) && ((memspace_out == memory::DEVICE)))
        {
          control = 3;
        }

        if (memspace_out == memory::HOST)
        {
          assert(((h_row_data_ == nullptr) == (h_col_data_ == nullptr)) && "In CsrMatrix::copyDataFrom one of host row or column data is null!\n");

          if ((h_row_data_ == nullptr) && (h_col_data_ == nullptr))
          {
            this->h_row_data_          = new IdxT[static_cast<size_t>(n_ + 1)];
            this->h_col_data_          = new IdxT[static_cast<size_t>(nnz_current)];
            owns_cpu_sparsity_pattern_ = true;
          }
          if (h_val_data_ == nullptr)
          {
            this->h_val_data_ = new RealT[static_cast<size_t>(nnz_current)];
            owns_cpu_values_  = true;
          }
        }

        if (memspace_out == memory::DEVICE)
        {
          assert(((d_row_data_ == nullptr) == (d_col_data_ == nullptr)) && "In CsrMatrix::copyDataFrom one of device row or column data is null!\n");

          if ((d_row_data_ == nullptr) && (d_col_data_ == nullptr))
          {
            mem_.allocateArrayOnDevice(&d_row_data_, n_ + 1);
            mem_.allocateArrayOnDevice(&d_col_data_, nnz_current);
            owns_gpu_values_ = true;
          }
          if (d_val_data_ == nullptr)
          {
            mem_.allocateArrayOnDevice(&d_val_data_, nnz_current);
            owns_gpu_sparsity_pattern_ = true;
          }
        }

        switch (control)
        {
        case 0: // cpu->cpu
          mem_.copyArrayHostToHost(h_row_data_, row_data, n_ + 1);
          mem_.copyArrayHostToHost(h_col_data_, col_data, nnz_current);
          mem_.copyArrayHostToHost(h_val_data_, val_data, nnz_current);
          h_data_updated_ = true;
          break;
        case 2: // gpu->cpu
          mem_.copyArrayDeviceToHost(h_row_data_, row_data, n_ + 1);
          mem_.copyArrayDeviceToHost(h_col_data_, col_data, nnz_current);
          mem_.copyArrayDeviceToHost(h_val_data_, val_data, nnz_current);
          h_data_updated_ = true;
          break;
        case 1: // cpu->gpu
          mem_.copyArrayHostToDevice(d_row_data_, row_data, n_ + 1);
          mem_.copyArrayHostToDevice(d_col_data_, col_data, nnz_current);
          mem_.copyArrayHostToDevice(d_val_data_, val_data, nnz_current);
          d_data_updated_ = true;
          break;
        case 3: // gpu->gpu
          mem_.copyArrayDeviceToDevice(d_row_data_, row_data, n_ + 1);
          mem_.copyArrayDeviceToDevice(d_col_data_, col_data, nnz_current);
          mem_.copyArrayDeviceToDevice(d_val_data_, val_data, nnz_current);
          d_data_updated_ = true;
          break;
        default:
          return -1;
        }
        return 0;
      }

      int copyDataFrom(const IdxT*         row_data,
                       const IdxT*         col_data,
                       const RealT*        val_data,
                       IdxT                new_nnz,
                       memory::MemorySpace memspace_in,
                       memory::MemorySpace memspace_out)
      {
        destroyMatrixData(memspace_out);
        nnz_ = new_nnz;
        return copyDataFrom(row_data, col_data, val_data, memspace_in, memspace_out);
      }

      /**
       * @brief Build CSR from sorted/deduplicated COO entries.
       *
       * Copies the provided row pointers, column indices, and values
       * into the host arrays of this matrix. Allocates host data if needed.
       *
       * @param[in] row_data - Row pointers (size n_ + 1)
       * @param[in] col_data - Column indices (size nnz_)
       * @param[in] val_data - Values (size nnz_)
       */
      void addEntries(const IdxT* row_data, const IdxT* col_data, const RealT* val_data)
      {
        for (IdxT i = 0; i < n_ + 1; ++i)
        {
          h_row_data_[i] = row_data[i];
        }
        for (IdxT k = 0; k < nnz_; ++k)
        {
          h_col_data_[k] = col_data[k];
          h_val_data_[k] = val_data[k];
        }
        h_data_updated_ = true;
        d_data_updated_ = false;
      }

      int allocateMatrixData(memory::MemorySpace memspace)
      {
        IdxT nnz_current = nnz_;
        destroyMatrixData(memspace);

        if (memspace == memory::HOST)
        {
          this->h_row_data_ = new IdxT[static_cast<size_t>(n_ + 1)];
          std::fill(h_row_data_, h_row_data_ + n_ + 1, 0);
          this->h_col_data_ = new IdxT[static_cast<size_t>(nnz_current)];
          std::fill(h_col_data_, h_col_data_ + nnz_current, 0);
          this->h_val_data_ = new RealT[static_cast<size_t>(nnz_current)];
          std::fill(h_val_data_, h_val_data_ + nnz_current, 0.0);
          owns_cpu_sparsity_pattern_ = true;
          owns_cpu_values_           = true;
          return 0;
        }

        if (memspace == memory::DEVICE)
        {
          mem_.allocateArrayOnDevice(&d_row_data_, n_ + 1);
          mem_.allocateArrayOnDevice(&d_col_data_, nnz_current);
          mem_.allocateArrayOnDevice(&d_val_data_, nnz_current);
          owns_gpu_sparsity_pattern_ = true;
          owns_gpu_values_           = true;
          return 0;
        }
        return -1;
      }

      /**
       * @brief Set the pointers for matrix row, column, value data.
       *
       * Useful if interfacing with other codes - this function only assigns
       * pointers, but it does not allocate nor copy anything.
       *
       * @param[in] row_data - pointer to row data
       * @param[in] col_data - pointer to column data
       * @param[in] val_data - pointer to value data
       * @param[in] memspace - memory space (HOST or DEVICE) of incoming data
       *
       * @return 0 if successful, 1 if not.
       */
      int setDataPointers(IdxT*               row_data,
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
       * @brief destroy matrix data (HOST or DEVICE) if the matrix owns it.
       *
       * @param[in] memspace - memory space (HOST or DEVICE)
       *
       * @return 0 if successful, -1 if not.
       */
      int destroyMatrixData(memory::MemorySpace memspace)
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
       * @brief Prints matrix data.
       *
       * @param out - Output stream where the matrix data is printed
       * @param indexing_base - base index for printing (default 0)
       */
      void print(std::ostream& out = std::cout, IdxT indexing_base = 0)
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

      /**
       * @brief Sync data in memspace with the updated memory space.
       *
       * @param memspace - memory space to be synced up (HOST or DEVICE)
       * @return int - 0 if successful, error code otherwise
       *
       * @pre The memory space other than `memspace` must be up-to-date.
       */
      int syncData(memory::MemorySpace memspace)
      {
        using namespace memory;

        switch (memspace)
        {
        case HOST:
          assert(((h_row_data_ == nullptr) == (h_col_data_ == nullptr)) && "In CsrMatrix::syncData one of host row or column data is null!\n");

          if (h_data_updated_)
          {
            std::cerr << "CsrMatrix::syncData is trying to sync host, but host already up to date!\n";
            assert(!h_data_updated_);
            return 1;
          }
          if (!d_data_updated_)
          {
            std::cerr << "CsrMatrix::syncData is trying to sync host with device, but device is out of date!\n"
                      << "See CsrMatrix::syncData documentation\n.";
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
          assert(((d_row_data_ == nullptr) == (d_col_data_ == nullptr)) && "In CsrMatrix::syncData one of device row or column data is null!\n");

          if (d_data_updated_)
          {
            std::cerr << "CsrMatrix::syncData is trying to sync device, but device already up to date!\n";
            assert(!d_data_updated_);
            return 1;
          }
          if (!h_data_updated_)
          {
            std::cerr << "CsrMatrix::syncData is trying to sync device with host, but host is out of date!\n"
                      << "See CsrMatrix::syncData documentation\n.";
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
        }
      }

      /**
       * @brief Update matrix values by copying from new_vals.
       *
       * Allocates if needed. Sets ownership and update flags.
       *
       * @param[in] new_vals    - pointer to new values data
       * @param[in] memspace_in  - memory space of new_vals
       * @param[in] memspace_out - memory space of matrix values to update
       *
       * @return 0 if successful, -1 if not.
       */
      int copyValues(const RealT*        new_vals,
                     memory::MemorySpace memspace_in,
                     memory::MemorySpace memspace_out)
      {
        IdxT nnz_current = nnz_;
        setNotUpdated();
        int control = -1;
        if ((memspace_in == memory::HOST) && (memspace_out == memory::HOST))
        {
          control = 0;
        }
        if ((memspace_in == memory::HOST) && (memspace_out == memory::DEVICE))
        {
          control = 1;
        }
        if ((memspace_in == memory::DEVICE) && (memspace_out == memory::HOST))
        {
          control = 2;
        }
        if ((memspace_in == memory::DEVICE) && (memspace_out == memory::DEVICE))
        {
          control = 3;
        }

        if (memspace_out == memory::HOST)
        {
          if (h_val_data_ == nullptr)
          {
            this->h_val_data_ = new RealT[static_cast<size_t>(nnz_current)];
            owns_cpu_values_  = true;
          }
        }

        if (memspace_out == memory::DEVICE)
        {
          if (d_val_data_ == nullptr)
          {
            mem_.allocateArrayOnDevice(&d_val_data_, nnz_current);
            owns_gpu_values_ = true;
          }
        }

        switch (control)
        {
        case 0: // cpu->cpu
          mem_.copyArrayHostToHost(h_val_data_, new_vals, nnz_current);
          h_data_updated_ = true;
          break;
        case 2: // cuda->cpu
          mem_.copyArrayDeviceToHost(h_val_data_, new_vals, nnz_current);
          h_data_updated_ = true;
          break;
        case 1: // cpu->cuda
          mem_.copyArrayHostToDevice(d_val_data_, new_vals, nnz_current);
          d_data_updated_ = true;
          break;
        case 3: // cuda->cuda
          mem_.copyArrayDeviceToDevice(d_val_data_, new_vals, nnz_current);
          d_data_updated_ = true;
          break;
        default:
          return -1;
        }
        return 0;
      }

      /**
       * @brief Set values pointer (does not copy). Use caution.
       *
       * @param[in] new_vals - pointer to new values data
       * @param[in] memspace - memory space (HOST or DEVICE) of new_vals
       *
       * @return 0 if successful, -1 if not.
       */
      int setValuesPointer(RealT*              new_vals,
                           memory::MemorySpace memspace)
      {
        using namespace memory;
        setNotUpdated();

        switch (memspace)
        {
        case HOST:
          if (owns_cpu_values_ && h_val_data_)
          {
            std::cerr << "Trying to set matrix host values, but the values already set!\n";
            std::cerr << "Ignoring setValuesPointer function call ...\n";
            return 1;
          }
          h_val_data_      = new_vals;
          h_data_updated_  = true;
          owns_cpu_values_ = false;
          break;
        case DEVICE:
          if (owns_gpu_values_ && d_val_data_)
          {
            std::cerr << "Trying to set matrix device values, but the values already set!\n";
            std::cerr << "Ignoring setValuesPointer function call ...\n";
            return 1;
          }
          d_val_data_      = new_vals;
          d_data_updated_  = true;
          owns_gpu_values_ = false;
          break;
        default:
          return -1;
        }
        return 0;
      }

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

      void setNotUpdated()
      {
        h_data_updated_ = false;
        d_data_updated_ = false;
      }

      // Data ownership flags
      bool owns_cpu_sparsity_pattern_{false}; ///< for row/col data
      bool owns_cpu_values_{false};           ///< for nonzero values

      bool owns_gpu_sparsity_pattern_{false}; ///< for row/col data
      bool owns_gpu_values_{false};           ///< for nonzero values

      MemoryHandler mem_; ///< Device memory manager object
    };
  } // namespace LinearAlgebra
} // namespace GridKit
