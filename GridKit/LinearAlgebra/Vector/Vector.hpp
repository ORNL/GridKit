#pragma once

#include <cassert>
#include <string>

#include <GridKit/MemoryUtilities/MemoryUtils.hpp>

namespace GridKit
{
  namespace LinearAlgebra
  {
    namespace detail
    {
      // NOTE: this namespace exists to reduce template instantiation + decouple
      //       code where it isn't necessary

      /**
       * @brief log failure for host data.
       */
      [[gnu::cold, gnu::noinline]]
      void logHostUnsyncFailure();

      /**
       * @brief log failure for device data.
       */
      [[gnu::cold, gnu::noinline]]
      void logDeviceUnsyncFailure();

      /**
       * @brief log failure for bounds check.
       */
      template <typename IdxT>
      [[gnu::cold, gnu::noinline]]
      void logBoundsCheckFailure(IdxT j, IdxT k);

      /**
       * @brief log failure for @ref Vector::setDataUpdated bounds check.
       */
      template <typename IdxT>
      [[gnu::cold, gnu::noinline]]
      void logUpdatedBoundsCheckFailure(IdxT j, IdxT k);
    } // namespace detail

    /**
     * @brief This class implements vectors (dense arrays) and multivectors and
     * some basic utilities (get size, allocate, set data, get data, etc).
     *
     *
     * What you need to know:
     *  - Multivectors are stored in one array, organized column-wise. A vector
     *    is a multivector of size 1.
     *  - There is a mirroring memory approach: the class has DEVICE and HOST
     *    data pointers. If needed,  only one (or none) can be used or allocated.
     *    Unless triggered directly or by other function, the data is NOT
     *    automatically updated between HOST and DEVICE.
     *  - Constructor DOES NOT allocate memory. This has to be done separately.
     *  - You can get (and set) "raw" data easily, if needed.
     *  - There is memory ownership utility - vector can own memory (separate
     *    flags for HOST and DEVICE) or not, depending on how it is used.
     *
     * @author Kasia Swirydowicz <kasia.swirydowicz@pnnl.gov>
     * @author Slaven Peles <peless@ornl.gov>
     */
    template <typename ScalarT, typename IdxT>
    class Vector
    {
    public:
      Vector()
        : Vector(0)
      {
      }

      Vector(IdxT n);
      Vector(IdxT n, IdxT k);
      ~Vector();

      Vector(const Vector&)            = delete;
      Vector(Vector&&)                 = delete;
      Vector& operator=(const Vector&) = delete;
      Vector& operator=(Vector&&)      = delete;

      int copyFromExternal(const ScalarT*      source,
                           memory::MemorySpace memspaceIn  = memory::HOST,
                           memory::MemorySpace memspaceOut = memory::HOST);
      int copyFromExternal(const Vector&       source,
                           memory::MemorySpace memspaceIn  = memory::HOST,
                           memory::MemorySpace memspaceOut = memory::HOST);

      // TODO: in the future, we should consider making the `getData` family of
      //       methods accept the memory space through a template parameter
      //       instead to guarantee the optimization implemented here through
      //       constant propagation

      /**
       * @brief get a pointer to HOST or DEVICE vector data.
       *
       * @param[in] memspace  - Memory space of the pointer (HOST or DEVICE)
       *
       * @return pointer to the vector data (HOST or DEVICE). In case of multivectors,
       * vectors are stored column-wise.
       *
       * @note This function gives you access to the pointer, not to a copy.
       * If you change the values using the pointer, the vector values will
       * change too. Make sure to use setDataUpdated function to set the update
       * flags correctly after changing the values.
       */
      [[gnu::always_inline]]
      inline ScalarT* getData(memory::MemorySpace memspace = memory::HOST)
      {
        using memory::DEVICE;
        using memory::HOST;

        switch (memspace)
        {
        case HOST:
          if (!getHostUpdated(0)) [[unlikely]]
          {
            detail::logHostUnsyncFailure();
            return nullptr;
          }
          return h_data_;
        case DEVICE:
          if (!getDeviceUpdated(0)) [[unlikely]]
          {
            detail::logDeviceUnsyncFailure();
            return nullptr;
          }
          return d_data_;
        default:
          return nullptr;
        }
      }

      /**
       * @brief Get a pointer to HOST or DEVICE data of a vector in a multivector.
       *
       * @param[in] j         - Index of a vector in multivector
       * @param[in] memspace  - Memory space of the pointer (HOST or DEVICE)
       *
       * @return Pointer to the _j_th vector data (HOST or DEVICE).
       *
       * @pre `j` < `k_`, i.e., `j` is smaller than the number of vectors.
       *
       * @note This function gives you access to the pointer, not to a copy.
       * If you change the values using the pointer, the vector values will
       * change too. Call setDataUpdated() to update the staleness flags.
       */
      [[gnu::always_inline]]
      ScalarT* getData(IdxT j, memory::MemorySpace memspace = memory::HOST)
      {
        using memory::DEVICE;
        using memory::HOST;

        if (k_ <= j) [[unlikely]]
        {
          detail::logBoundsCheckFailure<IdxT>(j, k_);
          return nullptr;
        }

        switch (memspace)
        {
        case HOST:
          if (!getHostUpdated(j)) [[unlikely]]
          {
            detail::logHostUnsyncFailure();
            return nullptr;
          }
          return &h_data_[j * n_size_];
        case DEVICE:
          if (!getDeviceUpdated(j)) [[unlikely]]
          {
            detail::logDeviceUnsyncFailure();
            return nullptr;
          }
          return &d_data_[j * n_size_];
        default:
          return nullptr;
        }
      }

      /**
       * @brief get a pointer to HOST or DEVICE vector data.
       *
       * @param[in] memspace  - Memory space of the pointer (HOST or DEVICE)
       *
       * @return pointer to the vector data (HOST or DEVICE). In case of multivectors,
       * vectors are stored column-wise.
       */
      [[gnu::always_inline]]
      const ScalarT* getData(memory::MemorySpace memspace = memory::HOST) const
      {
        using memory::DEVICE;
        using memory::HOST;

        switch (memspace)
        {
        case HOST:
          if (!getHostUpdated(0)) [[unlikely]]
          {
            detail::logHostUnsyncFailure();
            return nullptr;
          }
          return h_data_;
        case DEVICE:
          if (!getDeviceUpdated(0)) [[unlikely]]
          {
            detail::logDeviceUnsyncFailure();
            return nullptr;
          }
          return d_data_;
        default:
          return nullptr;
        }
      }

      /**
       * @brief Get a const pointer to HOST or DEVICE data of a vector in a multivector.
       *
       * @param[in] j         - Index of a vector in multivector
       * @param[in] memspace  - Memory space of the pointer (HOST or DEVICE)
       *
       * @return Const pointer to the _j_th vector data (HOST or DEVICE).
       *
       * @pre `j` < `k_`, i.e., `j` is smaller than the number of vectors.
       */
      [[gnu::always_inline]]
      const ScalarT* getData(IdxT j, memory::MemorySpace memspace = memory::HOST) const
      {
        using memory::DEVICE;
        using memory::HOST;

        if (k_ <= j) [[unlikely]]
        {
          detail::logBoundsCheckFailure<IdxT>(j, k_);
          return nullptr;
        }

        switch (memspace)
        {
        case HOST:
          if (!getHostUpdated(j)) [[unlikely]]
          {
            detail::logHostUnsyncFailure();
            return nullptr;
          }
          return &h_data_[j * n_size_];
        case DEVICE:
          if (!getDeviceUpdated(j)) [[unlikely]]
          {
            detail::logDeviceUnsyncFailure();
            return nullptr;
          }
          return &d_data_[j * n_size_];
        default:
          return nullptr;
        }
      }

      IdxT getCapacity() const;
      IdxT getSize() const;
      IdxT getNumVectors() const;

      /**
       * @brief Set the flag to indicate that the data (HOST or DEVICE) has been
       * updated.
       *
       * Use this function if you update vector elements by accessing the raw data
       * pointer.
       *
       * @param[in] memspace - Memory space (HOST or DEVICE)
       *
       * @warning This is an expert level method. Use only if you know what
       * you are doing.
       */
      [[gnu::always_inline]]
      int setDataUpdated(memory::MemorySpace memspace = memory::HOST)
      {
        using namespace memory;
        switch (memspace)
        {
        case HOST:
          setHostUpdated(true);
          setDeviceUpdated(false);
          break;
        case DEVICE:
          setHostUpdated(false);
          setDeviceUpdated(true);
          break;
        }
        return 0;
      }

      /**
       * @brief Set the flag to indicate that the data (HOST or DEVICE) for
       * vector `j` in the multivector has been updated.
       *
       * Use this function if you update vector elements by accessing the raw data
       * pointer.
       *
       * @param[in] memspace - Memory space (HOST or DEVICE)
       *
       * @warning This is an expert level method. Use only if you know what
       * you are doing.
       */
      [[gnu::always_inline]]
      int setDataUpdated(IdxT j, memory::MemorySpace memspace = memory::HOST)
      {
        using namespace memory;

        if (k_ <= j) [[unlikely]]
        {
          detail::logUpdatedBoundsCheckFailure(j, k_);
          return 1;
        }

        switch (memspace)
        {
        case HOST:
          setHostUpdated(j, true);
          setDeviceUpdated(j, false);
          break;
        case DEVICE:
          setDeviceUpdated(j, true);
          setHostUpdated(j, false);
          break;
        }
        return 0;
      }

      int setData(ScalarT* data, memory::MemorySpace memspace = memory::HOST);
      int setData(ScalarT* data, IdxT size, memory::MemorySpace memspace = memory::HOST);
      int aliasOf(Vector&             parent,
                  IdxT                offset,
                  IdxT                size,
                  memory::MemorySpace memspace = memory::HOST);
      int allocate(memory::MemorySpace memspace = memory::HOST);
      int setToZero(memory::MemorySpace memspace = memory::HOST);
      int setToZero(IdxT i, memory::MemorySpace memspace = memory::HOST);
      int setToConst(ScalarT C, memory::MemorySpace memspace = memory::HOST);
      int setToConst(IdxT i, ScalarT C, memory::MemorySpace memspace = memory::HOST);
      int syncData(memory::MemorySpace memspaceOut = memory::HOST);
      int syncData(IdxT j, memory::MemorySpace memspaceOut = memory::HOST);
      int resize(IdxT new_n_current);
      int copyToExternal(ScalarT*            dest,
                         IdxT                i,
                         memory::MemorySpace memspaceSrc = memory::HOST,
                         memory::MemorySpace memspaceDst = memory::HOST);
      int copyToExternal(ScalarT*            dest,
                         memory::MemorySpace memspaceSrc = memory::HOST,
                         memory::MemorySpace memspaceDst = memory::HOST);

    private:
      [[gnu::always_inline]]
      void setHostUpdated(bool is_updated)
      {
        if (k_ <= 1 && owns_update_flags_) [[likely]]
        {
          cpu_updated_single_ = is_updated;
          return;
        }

        std::fill(cpu_updated_, cpu_updated_ + k_, is_updated);
      }

      [[gnu::always_inline]]
      void setHostUpdated(IdxT j, bool is_updated)
      {
        if (k_ <= 1 && owns_update_flags_)
        {
          assert(j == 0);
          cpu_updated_single_ = is_updated;
          return;
        }

        cpu_updated_[j] = is_updated;
      }

      [[gnu::always_inline]]
      void setDeviceUpdated(bool is_updated)
      {
        if (k_ <= 1 && owns_update_flags_) [[likely]]
        {
          gpu_updated_single_ = is_updated;
          return;
        }

        std::fill(gpu_updated_, gpu_updated_ + k_, is_updated);
      }

      [[gnu::always_inline]]
      void setDeviceUpdated(IdxT j, bool is_updated)
      {
        if (k_ <= 1 && owns_update_flags_)
        {
          assert(j == 0);
          gpu_updated_single_ = is_updated;
          return;
        }

        gpu_updated_[j] = is_updated;
      }

      [[gnu::always_inline]]
      bool getHostUpdated(IdxT j) const
      {
        if (k_ <= 1 && owns_update_flags_) [[likely]]
        {
          assert(j == 0);
          return cpu_updated_single_;
        }

        return cpu_updated_[j];
      }

      [[gnu::always_inline]]
      bool getDeviceUpdated(IdxT j) const
      {
        if (k_ <= 1 && owns_update_flags_) [[likely]]
        {
          assert(j == 0);
          return gpu_updated_single_;
        }

        return gpu_updated_[j];
      }

      IdxT     n_capacity_{0};   ///< vector capacity
      IdxT     k_{0};            ///< number of vectors in multivector
      IdxT     n_size_{0};       ///< actual size of the vector
      ScalarT* d_data_{nullptr}; ///< DEVICE data array
      ScalarT* h_data_{nullptr}; ///< HOST data array

      /// device data flags
      union
      {
        /// when k > 1, store the flags on the heap
        bool* gpu_updated_{nullptr};

        /// when k == 1, store the flag in-place for locality
        bool gpu_updated_single_;
      };

      /// host data flags
      union
      {
        /// when k > 1, store the flags on the heap
        bool* cpu_updated_{nullptr};

        /// when k == 1, store the flag in-place for locality
        bool cpu_updated_single_;
      };

      bool owns_gpu_data_{true};     ///< data ownership flag for DEVICE data
      bool owns_cpu_data_{true};     ///< data ownership flag for HOST data
      bool owns_update_flags_{true}; ///< false when the update flags belong to an aliased parent

      MemoryManager mem_; ///< Device memory manager object
    };
  } // namespace LinearAlgebra
} // namespace GridKit
