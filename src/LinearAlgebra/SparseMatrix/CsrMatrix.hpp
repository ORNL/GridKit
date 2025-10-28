#pragma once

#include <cstdint>
#include <cstring>
#include <iostream>

namespace GridKit
{
  namespace LinearAlgebra
  {
    namespace memory
    {
      enum MemorySpace
      {
        HOST = 0,
        DEVICE
      };

      /**
       * @brief Class containing dummy functions when there is no GPU support.
       *
       * @author Slaven Peles <peless@ornl.gov>
       */
      struct Cpu
      {
        /**
         * @brief Dummy function to stand in when GPU support is not enabled.
         */
        static void deviceSynchronize()
        {
          // Nothing to synchronize
        }

        /**
         * @brief Dummy function to stand in when GPU support is not enabled.
         *
         * @return Allways return success!
         */
        static int getLastDeviceError()
        {
          // not on device, nothing to get
          return 0;
        }

        /**
         * @brief Dummy function to notify us something is wrong.
         *
         * This will be called only if GPU device support is not built, so
         * trying to access a device should indicate a bug in the code.
         *
         * @return Allways return failure!
         */
        static int deleteOnDevice(void* /* v */)
        {
          std::cerr << "Trying to delete on a GPU device, but GPU support not available.\n";
          return -1;
        }

        /**
         * @brief Dummy function to notify us something is wrong.
         *
         * This will be called only if GPU device support is not built, so
         * trying to access a device should indicate a bug in the code.
         *
         * @return Allways return failure!
         */
        template <typename I, typename T>
        static int allocateArrayOnDevice(T** /* v */, I /* n */)
        {
          std::cerr << "Trying to allocate on a GPU device, but GPU support not available.\n";
          return -1;
        }

        /**
         * @brief Dummy function to notify us something is wrong.
         *
         * This will be called only if GPU device support is not built, so
         * trying to access a device should indicate a bug in the code.
         *
         * @return Allways return failure!
         */
        template <typename I, typename T>
        static int allocateBufferOnDevice(T** /* v */, I /* n */)
        {
          std::cerr << "Trying to allocate on a GPU device, but GPU support not available.\n";
          return -1;
        }

        /**
         * @brief Dummy function to notify us something is wrong.
         *
         * This will be called only if GPU device support is not built, so
         * trying to access a device should indicate a bug in the code.
         *
         * @return Allways return failure!
         */
        template <typename I, typename T>
        static int setZeroArrayOnDevice(T* /* v */, I /* n */)
        {
          std::cerr << "Trying to initialize array on a GPU device, but GPU support not available.\n";
          return -1;
        }

        /**
         * @brief Dummy function to notify us something is wrong.
         *
         * This will be called only if GPU device support is not built, so
         * trying to access a device should indicate a bug in the code.
         *
         * @return Allways return failure!
         */
        template <typename I, typename T>
        static int setArrayToConstOnDevice(T* /* v */, T /* c */, I /* n */)
        {
          std::cerr << "Trying to initialize array on a GPU device, but GPU support not available.\n";
          return -1;
        }

        /**
         * @brief Dummy function to notify us something is wrong.
         *
         * This will be called only if GPU device support is not built, so
         * trying to access a device should indicate a bug in the code.
         *
         * @return Allways return failure!
         */
        template <typename I, typename T>
        static int copyArrayDeviceToHost(T* /* dst */, const T* /* src */, I /* n */)
        {
          std::cerr << "Trying to copy from a GPU device, but GPU support not available.\n";
          return -1;
        }

        /**
         * @brief Dummy function to notify us something is wrong.
         *
         * This will be called only if GPU device support is not built, so
         * trying to access a device should indicate a bug in the code.
         *
         * @return Allways return failure!
         */
        template <typename I, typename T>
        static int copyArrayDeviceToDevice(T* /* dst */, const T* /* src */, I /* n */)
        {
          std::cerr << "Trying to copy to a GPU device, but GPU support not available.\n";
          return -1;
        }

        template <typename I, typename T>
        static int copyArrayHostToDevice(T* /* dst */, const T* /* src */, I /* n */)
        {
          std::cerr << "Trying to copy to a GPU device, but GPU support not available.\n";
          return -1;
        }

      }; // struct Cuda
    }; // namespace memory

    /**
     * @class MemoryUtils
     *
     * @brief Provides basic memory allocation, free and copy functions.
     *
     * This class provedes abstractions for memory management functiosn for
     * different GPU programming models.
     *
     * @tparam Policy - Memory management policy (vendor specific)
     *
     * @author Slaven Peles <peless@ornl.gov>
     */
    template <class Policy>
    class MemoryUtils
    {
    public:
      MemoryUtils()  = default;
      ~MemoryUtils() = default;

      void deviceSynchronize();
      int  getLastDeviceError();
      int  deleteOnDevice(void* v);

      template <typename I, typename T>
      int allocateArrayOnDevice(T** v, I n);

      template <typename I, typename T>
      int allocateBufferOnDevice(T** v, I n);

      template <typename I, typename T>
      int setZeroArrayOnDevice(T* v, I n);

      template <typename I, typename T>
      int setArrayToConstOnDevice(T* v, T c, I n);

      template <typename I, typename T>
      int copyArrayDeviceToHost(T* dst, const T* src, I n);

      template <typename I, typename T>
      int copyArrayDeviceToDevice(T* dst, const T* src, I n);

      template <typename I, typename T>
      int copyArrayHostToDevice(T* dst, const T* src, I n);

      ///
      /// Methods implemented here are always needed
      ///

      template <typename I, typename T>
      int allocateArrayOnHost(T** v, I n)
      {
        std::size_t arraysize = static_cast<std::size_t>(n) * sizeof(T);
        *v                    = new T[arraysize];
        return *v == nullptr ? 1 : 0;
      }

      template <typename T>
      int deleteOnHost(T* v)
      {
        delete[] v;
        v = nullptr;
        return 0;
      }

      template <typename I, typename T>
      int copyArrayHostToHost(T* dst, const T* src, I n)
      {
        std::size_t arraysize = static_cast<std::size_t>(n) * sizeof(T);
        memcpy(dst, src, arraysize);
        return 0;
      }

      template <typename I, typename T>
      int setZeroArrayOnHost(T* v, I n)
      {
        std::size_t arraysize = static_cast<std::size_t>(n) * sizeof(T);
        memset(v, 0, arraysize);
        return 0;
      }

      template <typename I, typename T>
      int setArrayToConstOnHost(T* v, T c, I n)
      {
        for (I i = 0; i < n; ++i)
        {
          v[i] = c;
        }
        return 0;
      }
    };

    template <class Policy>
    void MemoryUtils<Policy>::deviceSynchronize()
    {
      Policy::deviceSynchronize();
    }

    template <class Policy>
    int MemoryUtils<Policy>::getLastDeviceError()
    {
      return Policy::getLastDeviceError();
    }

    template <class Policy>
    int MemoryUtils<Policy>::deleteOnDevice(void* v)
    {
      return Policy::deleteOnDevice(v);
    }

    template <class Policy>
    template <typename I, typename T>
    int MemoryUtils<Policy>::allocateArrayOnDevice(T** v, I n)
    {
      return Policy::template allocateArrayOnDevice<I, T>(v, n);
    }

    template <class Policy>
    template <typename I, typename T>
    int MemoryUtils<Policy>::allocateBufferOnDevice(T** v, I n)
    {
      return Policy::template allocateBufferOnDevice<I, T>(v, n);
    }

    template <class Policy>
    template <typename I, typename T>
    int MemoryUtils<Policy>::setZeroArrayOnDevice(T* v, I n)
    {
      return Policy::template setZeroArrayOnDevice<I, T>(v, n);
    }

    template <class Policy>
    template <typename I, typename T>
    int MemoryUtils<Policy>::setArrayToConstOnDevice(T* v, T c, I n)
    {
      return Policy::template setArrayToConstOnDevice<I, T>(v, c, n);
    }

    template <class Policy>
    template <typename I, typename T>
    int MemoryUtils<Policy>::copyArrayDeviceToHost(T* dst, const T* src, I n)
    {
      return Policy::template copyArrayDeviceToHost<I, T>(dst, src, n);
    }

    template <class Policy>
    template <typename I, typename T>
    int MemoryUtils<Policy>::copyArrayDeviceToDevice(T* dst, const T* src, I n)
    {
      return Policy::template copyArrayDeviceToDevice<I, T>(dst, src, n);
    }

    template <class Policy>
    template <typename I, typename T>
    int MemoryUtils<Policy>::copyArrayHostToDevice(T* dst, const T* src, I n)
    {
      return Policy::template copyArrayHostToDevice<I, T>(dst, src, n);
    }

    using real_type     = double;
    using index_type    = std::int32_t;
    using MemoryHandler = MemoryUtils<memory::Cpu>;

    // TODO: Port from ReSolve
    // #ifdef RESOLVE_USE_GPU

    // // Check if GPU support is enabled in Re::Solve and set appropriate device memory manager.
    // #if defined RESOLVE_USE_CUDA
    // #include <resolve/cuda/CudaMemory.hpp>
    //     using MemoryHandler = MemoryUtils<memory::Cuda>;
    // #elif defined RESOLVE_USE_HIP
    // #include <resolve/hip/HipMemory.hpp>
    //     using MemoryHandler = MemoryUtils<memory::Hip>;
    // #else
    // #error Unrecognized device, probably bug in CMake configuration
    // #endif

    // #else

    // If no GPU support is present, set device memory manager to a dummy object.
    using MemoryHandler = MemoryUtils<memory::Cpu>;

    // #endif

    class Sparse
    {
    public:
      /// Supported sparse matrix formats
      enum SparseFormat
      {
        NONE,
        TRIPLET,
        COMPRESSED_SPARSE_ROW,
        COMPRESSED_SPARSE_COLUMN
      };

    public:
      // basic constructor
      Sparse();
      Sparse(index_type n, index_type m, index_type nnz);
      Sparse(index_type n,
             index_type m,
             index_type nnz,
             bool       symmetric,
             bool       expanded);
      virtual ~Sparse();

      // accessors
      index_type   getNumRows();
      index_type   getNumColumns();
      index_type   getNnz();
      SparseFormat getSparseFormat() const;

      bool symmetric();
      bool expanded();
      void setSymmetric(bool symmetric);
      void setExpanded(bool expanded);
      void setNnz(index_type nnz_new); // for resetting when removing duplicates
      int  setUpdated(memory::MemorySpace what);

      virtual index_type* getRowData(memory::MemorySpace memspace) = 0;
      virtual index_type* getColData(memory::MemorySpace memspace) = 0;
      virtual real_type*  getValues(memory::MemorySpace memspace)  = 0;

      virtual int copyDataFrom(const index_type*   row_data,
                               const index_type*   col_data,
                               const real_type*    val_data,
                               memory::MemorySpace memspaceIn,
                               memory::MemorySpace memspaceOut) = 0;
      virtual int copyDataFrom(const index_type*   row_data,
                               const index_type*   col_data,
                               const real_type*    val_data,
                               index_type          new_nnz,
                               memory::MemorySpace memspaceIn,
                               memory::MemorySpace memspaceOut) = 0;

      virtual int allocateMatrixData(memory::MemorySpace memspace) = 0;
      int         setDataPointers(index_type*         row_data,
                                  index_type*         col_data,
                                  real_type*          val_data,
                                  memory::MemorySpace memspace);

      int destroyMatrixData(memory::MemorySpace memspace);

      virtual void print(std::ostream& file_out, index_type indexing_base) = 0;

      virtual int syncData(memory::MemorySpace memspaceOut) = 0;

      // update Values just updates values; it allocates if necessary.
      // values have the same dimensions between different formats
      virtual int copyValues(const real_type*    new_vals,
                             memory::MemorySpace memspaceIn,
                             memory::MemorySpace memspaceOut);

      // set new values just sets the pointer, use caution.
      virtual int setValuesPointer(real_type*          new_vals,
                                   memory::MemorySpace memspace);

    protected:
      SparseFormat sparse_format_{NONE}; ///< Matrix format
      index_type   n_{0};                ///< number of rows
      index_type   m_{0};                ///< number of columns
      index_type   nnz_{0};              ///< number of non-zeros

      bool is_symmetric_{false}; ///< symmetry flag
      bool is_expanded_{false};  ///< "expanded" flag

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

    class CsrMatrix : public Sparse
    {
    public:
      CsrMatrix();

      CsrMatrix(index_type n, index_type m, index_type nnz);

      CsrMatrix(index_type n,
                index_type m,
                index_type nnz,
                bool       symmetric,
                bool       expanded);

      CsrMatrix(index_type          n,
                index_type          m,
                index_type          nnz,
                bool                symmetric,
                bool                expanded,
                index_type**        rows,
                index_type**        cols,
                real_type**         vals,
                memory::MemorySpace memspaceSrc,
                memory::MemorySpace memspaceDst);

      ~CsrMatrix();

      virtual index_type* getRowData(memory::MemorySpace memspace);
      virtual index_type* getColData(memory::MemorySpace memspace);
      virtual real_type*  getValues(memory::MemorySpace memspace);

      virtual int copyDataFrom(const index_type*   row_data,
                               const index_type*   col_data,
                               const real_type*    val_data,
                               memory::MemorySpace memspaceIn,
                               memory::MemorySpace memspaceOut);
      virtual int copyDataFrom(const index_type*   row_data,
                               const index_type*   col_data,
                               const real_type*    val_data,
                               index_type          new_nnz,
                               memory::MemorySpace memspaceIn,
                               memory::MemorySpace memspaceOut);

      virtual int allocateMatrixData(memory::MemorySpace memspace);

      virtual void print(std::ostream& file_out = std::cout, index_type indexing_base = 0);

      virtual int syncData(memory::MemorySpace memspaceOut);
    };
  } // namespace LinearAlgebra
} // namespace GridKit