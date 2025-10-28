#include "CsrMatrix.hpp"

#include <cassert>
#include <iomanip>
#include <limits>

namespace GridKit
{
  namespace LinearAlgebra
  {
    /**
     * @brief empty constructor that does absolutely nothing
     */
    Sparse::Sparse()
    {
    }

    /**
     * @brief basic constructor. It DOES NOT allocate any memory!
     *
     * @param[in] n   - number of rows
     * @param[in] m   - number of columns
     * @param[in] nnz - number of non-zeros
     */
    Sparse::Sparse(index_type n,
                   index_type m,
                   index_type nnz)
      : n_{n},
        m_{m},
        nnz_{nnz}
    {
      this->is_symmetric_ = false;
      this->is_expanded_  = true; // default is a normal non-symmetric fully expanded matrix

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
     * @brief another basic constructor. It DOES NOT allocate any memory!
     *
     * @param[in] n         - number of rows
     * @param[in] m         - number of columns
     * @param[in] nnz       - number of non-zeros
     * @param[in] symmetric - true if symmetric, false if non-symmetric
     * @param[in] expanded  - true if expanded, false if not
     */
    Sparse::Sparse(index_type n,
                   index_type m,
                   index_type nnz,
                   bool       symmetric,
                   bool       expanded)
      : n_{n},
        m_{m},
        nnz_{nnz},
        is_symmetric_{symmetric},
        is_expanded_{expanded}
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
     * @brief destructor
     * */
    Sparse::~Sparse()
    {
      this->destroyMatrixData(memory::HOST);
      this->destroyMatrixData(memory::DEVICE);
    }

    /**
     * @brief set the matrix update flags to false (for both HOST and DEVICE).
     */
    void Sparse::setNotUpdated()
    {
      h_data_updated_ = false;
      d_data_updated_ = false;
    }

    /**
     * @brief get number of matrix rows
     *
     * @return number of matrix rows.
     */
    index_type Sparse::getNumRows()
    {
      return this->n_;
    }

    /**
     * @brief get number of matrix columns
     *
     * @return number of matrix columns.
     */
    index_type Sparse::getNumColumns()
    {
      return this->m_;
    }

    /**
     * @brief get number of non-zeros in the matrix.
     *
     * @return number of non-zeros.
     */
    index_type Sparse::getNnz()
    {
      return this->nnz_;
    }

    Sparse::SparseFormat Sparse::getSparseFormat() const
    {
      return sparse_format_;
    }

    /**
     * @brief check if matrix is symmetric.
     *
     * @return true if symmetric, false otherwise.
     */
    bool Sparse::symmetric()
    {
      return is_symmetric_;
    }

    /**
     * @brief check if (symmetric) matrix is expanded.
     *
     * @return true if expanded, false otherwise.
     */
    bool Sparse::expanded()
    {
      return is_expanded_;
    }

    /**
     * @brief Set matrix symmetry property
     *
     * @param[in] symmetric - true to set matrix to symmetric and false to set to non-symmetric
     */
    void Sparse::setSymmetric(bool symmetric)
    {
      this->is_symmetric_ = symmetric;
    }

    /**
     * @brief Set matrix "expanded" property
     *
     * @param[in] expanded - true to set matrix to expanded and false to set to not expanded
     */
    void Sparse::setExpanded(bool expanded)
    {
      this->is_expanded_ = expanded;
    }

    /**
     * @brief Set number of non-zeros.
     *
     * @param[in] nnz_new - new number of non-zeros
     */
    void Sparse::setNnz(index_type nnz_new)
    {
      this->nnz_ = nnz_new;
    }

    /**
     * @brief Tags `memspace` as updated.
     *
     * @param[in] memspace - memory space (HOST or DEVICE) to set to "updated"
     *
     * @return 0 if successful, -1 if not.
     *
     * The method sets the boolean flag indicating that the `memspace` is updated.
     * It automatically sets the other data mirror to non-updated. You would
     * use this function if you update matrix data by accessing its raw pointers.
     * In such case, the matrix has no way of knowing which data is most recent, so
     * you have to tell it.
     *
     * @warning This is an expert-level function. Use only if you know what you are
     * doing.
     *
     * @note If you want to set both DEVICE and HOST memory to the same value
     * use syncData function.
     */
    int Sparse::setUpdated(memory::MemorySpace memspace)
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
    int Sparse::setDataPointers(index_type*         row_data,
                                index_type*         col_data,
                                real_type*          val_data,
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
    int Sparse::destroyMatrixData(memory::MemorySpace memspace)
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
     * @brief updata matrix values using the _new_values_ provided either as HOST or as DEVICE array.
     *
     * This function will copy the data (not just assign a pointer) and allocate if needed.
     * It also sets ownership and update flags.
     *
     * @param[in] new_vals    - pointer to new values data (array of real numbers)
     * @param[in] memspaceIn  - memory space (HOST or DEVICE) of _new_vals_
     * @param[in] memspaceOut - memory space (HOST or DEVICE) of matrix values to be updated.
     *
     * @return 0 if successful, -1 if not.
     */
    int Sparse::copyValues(const real_type*    new_vals,
                           memory::MemorySpace memspaceIn,
                           memory::MemorySpace memspaceOut)
    {

      index_type nnz_current = nnz_;
      // four cases (for now)
      setNotUpdated();
      int control = -1;
      if ((memspaceIn == memory::HOST) && (memspaceOut == memory::HOST))
      {
        control = 0;
      }
      if ((memspaceIn == memory::HOST) && (memspaceOut == memory::DEVICE))
      {
        control = 1;
      }
      if ((memspaceIn == memory::DEVICE) && (memspaceOut == memory::HOST))
      {
        control = 2;
      }
      if ((memspaceIn == memory::DEVICE) && (memspaceOut == memory::DEVICE))
      {
        control = 3;
      }

      if (memspaceOut == memory::HOST)
      {
        // check if cpu data allocated
        if (h_val_data_ == nullptr)
        {
          this->h_val_data_ = new real_type[nnz_current];
          owns_cpu_values_  = true;
        }
      }

      if (memspaceOut == memory::DEVICE)
      {
        // check if cuda data allocated
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
     * @brief updata matrix values using the _new_values_ provided either as
     * HOST or as DEVICE array.
     *
     * This function only assigns a pointer, but does not copy. It sets update
     * flags.
     *
     * @param[in] new_vals - pointer to new values data (array of real numbers)
     * @param[in] memspace - memory space (HOST or DEVICE) of _new_vals_
     *
     * @return 0 if successful, -1 if not.
     */
    int Sparse::setValuesPointer(real_type*          new_vals,
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

    CsrMatrix::CsrMatrix()
    {
      sparse_format_ = COMPRESSED_SPARSE_ROW;
    }

    CsrMatrix::CsrMatrix(index_type n, index_type m, index_type nnz)
      : Sparse(n, m, nnz)
    {
      sparse_format_ = COMPRESSED_SPARSE_ROW;
    }

    CsrMatrix::CsrMatrix(index_type n,
                         index_type m,
                         index_type nnz,
                         bool       symmetric,
                         bool       expanded)
      : Sparse(n, m, nnz, symmetric, expanded)
    {
      sparse_format_ = COMPRESSED_SPARSE_ROW;
    }

    /**
     * @brief Hijacking constructor
     *
     * @param[in] n
     * @param[in] m
     * @param[in] nnz
     * @param[in] symmetric
     * @param[in] expanded
     * @param[in,out] rows
     * @param[in,out] cols
     * @param[in,out] vals
     * @param[in] memspaceSrc
     * @param[in] memspaceDst
     */
    CsrMatrix::CsrMatrix(index_type          n,
                         index_type          m,
                         index_type          nnz,
                         bool                symmetric,
                         bool                expanded,
                         index_type**        rows,
                         index_type**        cols,
                         real_type**         vals,
                         memory::MemorySpace memspaceSrc,
                         memory::MemorySpace memspaceDst)
      : Sparse(n, m, nnz, symmetric, expanded)
    {
      sparse_format_ = COMPRESSED_SPARSE_ROW;

      int control = -1;
      if ((memspaceSrc == memory::HOST) && (memspaceDst == memory::HOST))
      {
        control = 0;
      }
      if ((memspaceSrc == memory::HOST) && (memspaceDst == memory::DEVICE))
      {
        control = 1;
      }
      if ((memspaceSrc == memory::DEVICE) && (memspaceDst == memory::HOST))
      {
        control = 2;
      }
      if ((memspaceSrc == memory::DEVICE) && (memspaceDst == memory::DEVICE))
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
        syncData(memspaceDst);
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
        syncData(memspaceDst);

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
        std::cerr << "CsrMatrix constructor failed! "
                  << "Possible bug in memory spaces setting.\n";
        break;
      }
    }

    CsrMatrix::~CsrMatrix()
    {
    }

    index_type* CsrMatrix::getRowData(memory::MemorySpace memspace)
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

    index_type* CsrMatrix::getColData(memory::MemorySpace memspace)
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

    real_type* CsrMatrix::getValues(memory::MemorySpace memspace)
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

    int CsrMatrix::copyDataFrom(const index_type*   row_data,
                                const index_type*   col_data,
                                const real_type*    val_data,
                                memory::MemorySpace memspaceIn,
                                memory::MemorySpace memspaceOut)
    {
      // four cases (for now)
      index_type nnz_current = nnz_;
      setNotUpdated();
      int control = -1;
      if ((memspaceIn == memory::HOST) && (memspaceOut == memory::HOST))
      {
        control = 0;
      }
      if ((memspaceIn == memory::HOST) && ((memspaceOut == memory::DEVICE)))
      {
        control = 1;
      }
      if (((memspaceIn == memory::DEVICE)) && (memspaceOut == memory::HOST))
      {
        control = 2;
      }
      if (((memspaceIn == memory::DEVICE)) && ((memspaceOut == memory::DEVICE)))
      {
        control = 3;
      }

      if (memspaceOut == memory::HOST)
      {
        // check if cpu data allocated
        assert(((h_row_data_ == nullptr) == (h_col_data_ == nullptr)) && "In CsrMatrix::copyDataFrom one of host row or column data is null!\n");

        if ((h_row_data_ == nullptr) && (h_col_data_ == nullptr))
        {
          this->h_row_data_          = new index_type[n_ + 1];
          this->h_col_data_          = new index_type[nnz_current];
          owns_cpu_sparsity_pattern_ = true;
        }
        if (h_val_data_ == nullptr)
        {
          this->h_val_data_ = new real_type[nnz_current];
          owns_cpu_values_  = true;
        }
      }

      if (memspaceOut == memory::DEVICE)
      {
        // check if cuda data allocated
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

      // copy
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

    int CsrMatrix::copyDataFrom(const index_type*   row_data,
                                const index_type*   col_data,
                                const real_type*    val_data,
                                index_type          new_nnz,
                                memory::MemorySpace memspaceIn,
                                memory::MemorySpace memspaceOut)
    {
      destroyMatrixData(memspaceOut);
      nnz_ = new_nnz;
      return copyDataFrom(row_data, col_data, val_data, memspaceIn, memspaceOut);
    }

    int CsrMatrix::allocateMatrixData(memory::MemorySpace memspace)
    {
      index_type nnz_current = nnz_;
      destroyMatrixData(memspace); // just in case

      if (memspace == memory::HOST)
      {
        this->h_row_data_ = new index_type[n_ + 1];
        std::fill(h_row_data_, h_row_data_ + n_ + 1, 0);
        this->h_col_data_ = new index_type[nnz_current];
        std::fill(h_col_data_, h_col_data_ + nnz_current, 0);
        this->h_val_data_ = new real_type[nnz_current];
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
     * @brief Sync data in memspace with the updated memory space.
     *
     * @param memspace - memory space to be synced up (HOST or DEVICE)
     * @return int - 0 if successful, error code otherwise
     *
     * @pre The memory space other than `memspace` must be up-to-date. Otherwise,
     * this function will return an error.
     *
     * @see Sparse::setUpdated
     */
    int CsrMatrix::syncData(memory::MemorySpace memspace)
    {
      using namespace memory;

      switch (memspace)
      {
      case HOST:
        // check if we need to copy or not
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
          h_row_data_                = new index_type[n_ + 1];
          h_col_data_                = new index_type[nnz_];
          owns_cpu_sparsity_pattern_ = true;
        }
        if (h_val_data_ == nullptr)
        {
          h_val_data_      = new real_type[nnz_];
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
      } // switch
    }

    /**
     * @brief Prints matrix data.
     *
     * @param out - Output stream where the matrix data is printed
     */
    void CsrMatrix::print(std::ostream& out, index_type indexing_base)
    {
      out << std::scientific << std::setprecision(std::numeric_limits<real_type>::digits10);
      for (index_type i = 0; i < n_; ++i)
      {
        for (index_type j = h_row_data_[i]; j < h_row_data_[i + 1]; ++j)
        {
          out << i + indexing_base << " "
              << h_col_data_[j] + indexing_base << " "
              << h_val_data_[j] << "\n";
        }
      }
    }
  } // namespace LinearAlgebra
} // namespace GridKit
