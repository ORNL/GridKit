#include <cassert>
#include <cstring>

#include <GridKit/LinearAlgebra/Vector/Vector.hpp>
#include <GridKit/Utilities/Logger/Logger.hpp>

namespace GridKit
{
  namespace LinearAlgebra
  {

    using out = GridKit::Utilities::Logger;

    /**
     * @brief Single vector constructor.
     *
     * @param[in] n - Number of elements in the vector
     */
    template <typename ScalarT, typename IdxT>
    Vector<ScalarT, IdxT>::Vector(IdxT n)
      : n_capacity_(n),
        k_(1),
        n_size_(n),
        gpu_updated_(new bool[1]),
        cpu_updated_(new bool[1])
    {
      gpu_updated_[0] = false;
      cpu_updated_[0] = false;
    }

    /**
     * @brief Multivector constructor.
     *
     * @param[in] n - Number of elements in the vector
     * @param[in] k - Number of vectors in multivector
     */
    template <typename ScalarT, typename IdxT>
    Vector<ScalarT, IdxT>::Vector(IdxT n, IdxT k)
      : n_capacity_(n),
        k_(k),
        n_size_(n),
        gpu_updated_(new bool[k]),
        cpu_updated_(new bool[k])
    {
      setHostUpdated(false);
      setDeviceUpdated(false);
    }

    /**
     * @brief destructor.
     *
     */
    template <typename ScalarT, typename IdxT>
    Vector<ScalarT, IdxT>::~Vector()
    {
      if (owns_cpu_data_ && h_data_)
        mem_.deleteOnHost(h_data_);
      if (owns_gpu_data_ && d_data_)
        mem_.deleteOnDevice(d_data_);
      delete[] gpu_updated_;
      delete[] cpu_updated_;
    }

    /**
     * @brief Get capacity of a single vector.
     *
     * Vector memory is allocated to `n_capacity_*k_`. This is the maximum
     * number of elements that the (multi)vector can hold.
     *
     * @return `n_capacity_` the maximum number of elements in the vector.
     */
    template <typename ScalarT, typename IdxT>
    IdxT Vector<ScalarT, IdxT>::getCapacity() const
    {
      return n_capacity_;
    }

    /**
     * @brief get the number of elements in a single vector.
     *
     * For vectors with changing sizes, set the vector capacity to
     * the maximum expected size.
     *
     * @return `n_size_` number of elements currently in the vector.
     */
    template <typename ScalarT, typename IdxT>
    IdxT Vector<ScalarT, IdxT>::getSize() const
    {
      return n_size_;
    }

    /**
     * @brief Get the number of vectors in multivector.
     *
     * @return _k_, number of vectors in the multivector,
     * or 1 if the vector is not a multivector.
     */
    template <typename ScalarT, typename IdxT>
    IdxT Vector<ScalarT, IdxT>::getNumVectors() const
    {
      return k_;
    }

    /**
     * @brief Set the vector data pointer (HOST or DEVICE) to an external data.
     *
     * @param[in] data     - Pointer to data
     * @param[in] memspace - Memory space (HOST or DEVICE)
     *
     * @pre Vector data pointers must be null. If the vector data already
     * exists this function returns error message.
     *
     * @warning This function DOES NOT ALLOCATE any data, it only assigns the
     * pointer.
     *
     * @warning This is an expert level method. Use only if you know what
     * you are doing.
     */
    template <typename ScalarT, typename IdxT>
    int Vector<ScalarT, IdxT>::setData(ScalarT* data, memory::MemorySpace memspace)
    {
      using namespace memory;
      if (d_data_ || h_data_)
      {
        out::error() << "Trying to set vector data, but the data already exists!\n";
        out::error() << "Ignoring setData function call ...\n";
        return 1;
      }

      switch (memspace)
      {
      case HOST:
        h_data_ = data;
        setHostUpdated(true);
        setDeviceUpdated(false);
        owns_cpu_data_ = false;
        break;
      case DEVICE:
        d_data_ = data;
        setHostUpdated(false);
        setDeviceUpdated(true);
        owns_gpu_data_ = false;
        break;
      }
      return 0;
    }

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
    template <typename ScalarT, typename IdxT>
    int Vector<ScalarT, IdxT>::setDataUpdated(memory::MemorySpace memspace)
    {
      assert(cpu_updated_ && gpu_updated_ && "Update flags not allocated");

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
    template <typename ScalarT, typename IdxT>
    int Vector<ScalarT, IdxT>::setDataUpdated(IdxT j, memory::MemorySpace memspace)
    {
      assert(cpu_updated_ && gpu_updated_ && "Update flags not allocated");

      using namespace memory;
      switch (memspace)
      {
      case HOST:
        cpu_updated_[j] = true;
        gpu_updated_[j] = false;
        break;
      case DEVICE:
        gpu_updated_[j] = true;
        cpu_updated_[j] = false;
        break;
      }
      return 0;
    }

    /**
     * @brief Copy data from another vector.
     *
     * @param[in] v           - Vector, which data will be copied
     * @param[in] memspaceSrc  - Memory space of the data source (HOST or DEVICE)
     * @param[in] memspaceDst - Memory space the data will be copied to (HOST or DEVICE)
     *
     * @pre   size of _v_ is equal or larger than the current vector size.
     */
    template <typename ScalarT, typename IdxT>
    int Vector<ScalarT, IdxT>::copyFromExternal(const Vector& source, memory::MemorySpace memspaceSrc, memory::MemorySpace memspaceDst)
    {
      const ScalarT* source_data = source.getData(memspaceSrc);
      return copyFromExternal(source_data, memspaceSrc, memspaceDst);
    }

    /**
     * @brief Copy vector data from input array.
     *
     * This function allocates (if necessary) and copies the data.
     *
     * @param[in] data        - Data that is to be copied
     * @param[in] memspaceSrc - Memory space of the data source (HOST or DEVICE)
     * @param[in] memspaceDst - Memory space the data will be copied to (HOST or DEVICE)
     *
     * @return 0 if successful, 1 otherwise.
     */
    template <typename ScalarT, typename IdxT>
    int Vector<ScalarT, IdxT>::copyFromExternal(const ScalarT*      source,
                                                memory::MemorySpace memspaceSrc,
                                                memory::MemorySpace memspaceDst)
    {
      switch (memspaceDst)
      {
      case memory::HOST:
        if (h_data_ == nullptr)
        {
          out::error() << "Vector::copyFromExternal - destination memory space"
                       << " on host not allocated\n";
          return 1;
        }
        break;
      case memory::DEVICE:
        if (d_data_ == nullptr)
        {
          out::error() << "Vector::copyFromExternal - destination memory space"
                       << " on device not allocated\n";
          return 1;
        }
        break;
      }

      switch (memspaceSrc)
      {
      case memory::HOST:
        switch (memspaceDst)
        {
        case memory::HOST:
          mem_.copyArrayHostToHost(h_data_, source, n_size_ * k_);
          setHostUpdated(true);
          setDeviceUpdated(false);
          break;
        case memory::DEVICE:
          mem_.copyArrayHostToDevice(d_data_, source, n_size_ * k_);
          setHostUpdated(false);
          setDeviceUpdated(true);
          break;
        default:
          return 1;
        }
        break;
      case memory::DEVICE:
        switch (memspaceDst)
        {
        case memory::HOST:
          mem_.copyArrayDeviceToHost(h_data_, source, n_size_ * k_);
          setHostUpdated(true);
          setDeviceUpdated(false);
          break;
        case memory::DEVICE:
          mem_.copyArrayDeviceToDevice(d_data_, source, n_size_ * k_);
          setHostUpdated(false);
          setDeviceUpdated(true);
          break;
        default:
          return 1;
        }
        break;
      default:
        return 1;
      }
      return 0;
    }

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
    template <typename ScalarT, typename IdxT>
    ScalarT* Vector<ScalarT, IdxT>::getData(memory::MemorySpace memspace)
    {
      using memory::DEVICE;
      using memory::HOST;

      switch (memspace)
      {
      case HOST:
        if (cpu_updated_[0] == false)
        {
          out::error() << "Vector::getData - requesting data from host"
                       << " but host data is stale. Returning nullptr ...\n";
          return nullptr;
        }
        return h_data_;
      case DEVICE:
        if (gpu_updated_[0] == false)
        {
          out::error() << "Vector::getData - requesting data from device"
                       << " but device data is stale. Returning nullptr ...\n";
          return nullptr;
        }
        return d_data_;
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
    template <typename ScalarT, typename IdxT>
    const ScalarT* Vector<ScalarT, IdxT>::getData(memory::MemorySpace memspace) const
    {
      using memory::DEVICE;
      using memory::HOST;

      switch (memspace)
      {
      case HOST:
        if (cpu_updated_[0] == false)
        {
          out::error() << "Trying to get data on the host, but host data is out of date!\n"
                       << "Use syncData function to sync host data with the device data!\n";
          return nullptr;
        }
        return h_data_;
      case DEVICE:
        if (gpu_updated_[0] == false)
        {
          out::error() << "Trying to get data on the device, but device data is out of date!\n"
                       << "Use syncData function to sync device data with the host data!\n";
          return nullptr;
        }
        return d_data_;
      default:
        return nullptr;
      }
    }

    /**
     * @brief get a pointer to HOST or DEVICE data of a particular vector in a multivector.
     *
     * @param[in] j         - Index of a vector in multivector
     * @param[in] memspace  - Memory space of the pointer (HOST or DEVICE)
     *
     * @return pointer to the _i_th vector data (HOST or DEVICE) within a multivector.
     *
     * @pre `j` < `k_` i.e, `j` is smaller than the total number of vectors in multivector.
     *
     * @note This function gives you access to the pointer, not to a copy.
     * If you change the values using the pointer, the vector values will
     * change too. Make sure to use setDataUpdated function to set the update
     * flags correctly after changing the values.
     */
    template <typename ScalarT, typename IdxT>
    ScalarT* Vector<ScalarT, IdxT>::getData(IdxT j, memory::MemorySpace memspace)
    {
      using memory::DEVICE;
      using memory::HOST;

      if (k_ <= j)
      {
        out::error() << "Trying to get data for vector " << j << " in multivector"
                     << " but there are only " << k_ << " vectors!\n";
        return nullptr;
      }

      switch (memspace)
      {
      case HOST:
        if (cpu_updated_[j] == false)
        {
          out::error() << "Trying to get data for vector " << j << " on the host, "
                       << "but host data is out of date!\n"
                       << "Use syncData function to sync host data with the device data!\n";
          return nullptr;
        }
        return &h_data_[j * n_size_];
      case DEVICE:
        if (gpu_updated_[j] == false)
        {
          out::error() << "Trying to get data for vector " << j << " on the device, "
                       << "but device data is out of date!\n"
                       << "Use syncData function to sync device data with the host data!\n";
          return nullptr;
        }
        return &d_data_[j * n_size_];
      default:
        return nullptr;
      }
    }

    /**
     * @brief get a const pointer to HOST or DEVICE data of a particular
     * vector in a multivector.
     *
     * @param[in] j         - Index of a vector in multivector
     * @param[in] memspace  - Memory space of the pointer (HOST or DEVICE)
     *
     * @return pointer to the _i_th vector data (HOST or DEVICE) within a multivector.
     *
     * @pre `j` < `k_` i.e, `j` is smaller than the total number of vectors in multivector.
     *
     */
    template <typename ScalarT, typename IdxT>
    const ScalarT* Vector<ScalarT, IdxT>::getData(IdxT j, memory::MemorySpace memspace) const
    {
      using memory::DEVICE;
      using memory::HOST;

      if (k_ <= j)
      {
        out::error() << "Trying to get data for vector " << j << " in multivector"
                     << " but there are only " << k_ << " vectors!\n";
        return nullptr;
      }

      switch (memspace)
      {
      case HOST:
        if (cpu_updated_[j] == false)
        {
          out::error() << "Trying to get data for vector " << j << " on the host, "
                       << "but host data is out of date!\n"
                       << "Use syncData function to sync host data with the device data!\n";
          return nullptr;
        }
        return &h_data_[j * n_size_];
      case DEVICE:
        if (gpu_updated_[j] == false)
        {
          out::error() << "Trying to get data for vector " << j << " on the device, "
                       << "but device data is out of date!\n"
                       << "Use syncData function to sync device data with the host data!\n";
          return nullptr;
        }
        return &d_data_[j * n_size_];
      default:
        return nullptr;
      }
    }

    /**
     * @brief Sync out of date memory space with the updated one.
     *
     * syncData is the only function that can set data on both HOST and DEVICE
     * to the same values.
     *
     * @param[in] memspaceDst  - Memory space to sync
     *
     * @return 0 if successful, 1 otherwise.
     *
     * @warning This function can be called only when all vectors in a
     * multivector have the same update status. Otherwise, you need to sync
     * vectors in a multivector individually.
     *
     */
    template <typename ScalarT, typename IdxT>
    int Vector<ScalarT, IdxT>::syncData(memory::MemorySpace memspaceDst)
    {
      using namespace memory;

      bool all_cpu_updated = cpu_updated_[0];
      bool all_gpu_updated = gpu_updated_[0];

      // Verify that all vectors in multivector have the same update status.
      for (IdxT i = 1; i < k_; ++i)
      {
        if (gpu_updated_[i] != all_gpu_updated)
        {
          out::error() << "Trying to sync all multivector data on the device,"
                       << " but individual vectors were updated differently.\n"
                       << "Use syncData function for individual vectors instead!\n";
          return 1;
        }
        if (cpu_updated_[i] != all_cpu_updated)
        {
          out::error() << "Trying to sync all multivector data on the host,"
                       << " but individual vectors were updated differently.\n"
                       << "Use syncData function for individual vectors instead!\n";
          return 1;
        }
      }

      switch (memspaceDst)
      {
      case DEVICE: // cpu -> gpu
        if (gpu_updated_[0])
        {
          out::error() << "Trying to sync device, but device already up to date!\n";
          return 1;
        }
        if (!cpu_updated_[0])
        {
          out::error() << "Trying to sync device with host, but host is out of date!\n";
          return 1;
        }
        if (d_data_ == nullptr)
        {
          out::error() << "Vector::syncData - device data not allocated\n";
          return 1;
        }
        mem_.copyArrayHostToDevice(d_data_, h_data_, n_size_ * k_);
        setDeviceUpdated(true);
        break;
      case HOST: // gpu -> cpu
        if (cpu_updated_[0])
        {
          out::error() << "Trying to sync host, but host already up to date!\n";
          return 1;
        }
        if (!gpu_updated_[0])
        {
          out::error() << "Trying to sync host with device, but device is out of date!\n";
          return 1;
        }
        if (h_data_ == nullptr)
        {
          out::error() << "Vector::syncData - host data not allocated\n";
          return 1;
        }
        mem_.copyArrayDeviceToHost(h_data_, d_data_, n_size_ * k_);
        setHostUpdated(true);
        break;
      default:
        return 1;
      }
      return 0;
    }

    /**
     * @brief Sync out of date memory space with the updated one.
     *
     * syncData is the only function that can set data on both HOST and DEVICE
     * to the same values.
     *
     * @param[in] memspaceDst  - Memory space to sync
     *
     * @return 0 if successful, 1 otherwise.
     *
     * @warning This function can be called only when all vectors in a
     * multivector have the same update status. Otherwise, you need to sync
     * vectors in a multivector individually.
     *
     */
    template <typename ScalarT, typename IdxT>
    int Vector<ScalarT, IdxT>::syncData(IdxT j, memory::MemorySpace memspaceDst)
    {
      using namespace memory;

      switch (memspaceDst)
      {
      case DEVICE: // cpu->gpu
        if (gpu_updated_[j])
        {
          out::error() << "Trying to sync device, but device already up to date!\n";
          return 1;
        }
        if (!cpu_updated_[j])
        {
          out::error() << "Trying to sync device with host, but host is out of date!\n";
          return 1;
        }
        if (d_data_ == nullptr)
        {
          out::error() << "Vector::syncData - device data not allocated\n";
          return 1;
        }
        mem_.copyArrayHostToDevice(&d_data_[j * n_size_], &h_data_[j * n_size_], n_size_);
        gpu_updated_[j] = true;
        break;
      case HOST: // cuda->cpu
        if (cpu_updated_[j])
        {
          out::error() << "Trying to sync host, but host already up to date!\n";
          return 1;
        }
        if (!gpu_updated_[j])
        {
          out::error() << "Trying to sync host with device, but device is out of date!\n";
          return 1;
        }
        if (h_data_ == nullptr)
        {
          out::error() << "Vector::syncData - host data not allocated\n";
          return 1;
        }
        mem_.copyArrayDeviceToHost(&h_data_[j * n_size_], &d_data_[j * n_size_], n_size_);
        cpu_updated_[j] = true;
        break;
      default:
        return 1;
      }
      return 0;
    }

    /**
     * @brief Allocate vector data for HOST or DEVICE
     *
     * @param[in] memspace   - Memory space of the data to be allocated
     *
     */
    template <typename ScalarT, typename IdxT>
    int Vector<ScalarT, IdxT>::allocate(memory::MemorySpace memspace)
    {
      using namespace memory;
      switch (memspace)
      {
      case HOST:
        if (!owns_cpu_data_)
        {
          out::error() << "Vector::allocate - cannot allocate host data, vector does not own it\n";
          return 1;
        }
        mem_.deleteOnHost(h_data_);
        mem_.allocateArrayOnHost(&h_data_, n_capacity_ * k_);
        owns_cpu_data_ = true;
        setHostUpdated(false);
        break;
      case DEVICE:
        if (!owns_gpu_data_)
        {
          out::error() << "Vector::allocate - cannot allocate device data, vector does not own it\n";
          return 1;
        }
        mem_.deleteOnDevice(d_data_);
        mem_.allocateArrayOnDevice(&d_data_, n_capacity_ * k_);
        owns_gpu_data_ = true;
        setDeviceUpdated(false);
        break;
      }
      return 0;
    }

    /**
     * @brief set vector data to zero. In case of multivectors, entire multivector is set to zero.
     *
     * @param[in] memspace   - Memory space of the data to be set to 0 (HOST or DEVICE)
     *
     */
    template <typename ScalarT, typename IdxT>
    int Vector<ScalarT, IdxT>::setToZero(memory::MemorySpace memspace)
    {
      using namespace memory;
      switch (memspace)
      {
      case HOST:
        if (h_data_ == nullptr)
        {
          out::error() << "Vector::setToZero - host data not allocated\n";
          return 1;
        }
        mem_.setZeroArrayOnHost(h_data_, n_capacity_ * k_);
        setHostUpdated(true);
        setDeviceUpdated(false);
        break;
      case DEVICE:
        if (d_data_ == nullptr)
        {
          out::error() << "Vector::setToZero - device data not allocated\n";
          return 1;
        }
        mem_.setZeroArrayOnDevice(d_data_, n_capacity_ * k_);
        setHostUpdated(false);
        setDeviceUpdated(true);
        break;
      }
      return 0;
    }

    /**
     * @brief set the data of a single vector in a multivector to zero.
     *
     * @param[in] i          - Index of a vector in a multivector
     * @param[in] memspace   - Memory space of the data to be set to 0 (HOST or DEVICE)
     *
     * @pre   _i_ < _k_ i.e,, _i_ is smaller than the total number of vectors in multivector.
     */
    template <typename ScalarT, typename IdxT>
    int Vector<ScalarT, IdxT>::setToZero(IdxT j, memory::MemorySpace memspace)
    {
      using namespace memory;
      switch (memspace)
      {
      case HOST:
        if (h_data_ == nullptr)
        {
          out::error() << "Vector::setToZero - host data not allocated\n";
          return 1;
        }
        mem_.setZeroArrayOnHost(&h_data_[j * n_size_], n_size_);
        cpu_updated_[j] = true;
        gpu_updated_[j] = false;
        break;
      case DEVICE:
        if (d_data_ == nullptr)
        {
          out::error() << "Vector::setToZero - device data not allocated\n";
          return 1;
        }
        // TODO: We should not need to access raw data in this class
        mem_.setZeroArrayOnDevice(&d_data_[j * n_size_], n_size_);
        cpu_updated_[j] = false;
        gpu_updated_[j] = true;
        break;
      }
      return 0;
    }

    /**
     * @brief set vector data to a given constant.
     *
     * In case of multivectors, entire multivector is set to the constant.
     *
     * @param[in] C          - Constant (real number)
     * @param[in] memspace   - Memory space of the data to be set to 0 (HOST or DEVICE)
     *
     */
    template <typename ScalarT, typename IdxT>
    int Vector<ScalarT, IdxT>::setToConst(ScalarT C, memory::MemorySpace memspace)
    {
      using namespace memory;
      switch (memspace)
      {
      case HOST:
        if (h_data_ == nullptr)
        {
          out::error() << "Vector::setToConst - host data not allocated\n";
          return 1;
        }
        mem_.setArrayToConstOnHost(h_data_, C, n_size_ * k_);
        setHostUpdated(true);
        setDeviceUpdated(false);
        break;
      case DEVICE:
        if (d_data_ == nullptr)
        {
          out::error() << "Vector::setToConst - device data not allocated\n";
          return 1;
        }
        mem_.setArrayToConstOnDevice(d_data_, C, n_size_ * k_);
        setHostUpdated(false);
        setDeviceUpdated(true);
        break;
      }
      return 0;
    }

    /**
     * @brief set the data of a single vector in a multivector to a given constant.
     *
     * @param[in] j          - Index of a vector in a multivector
     * @param[in] C          - Constant (real number)
     * @param[in] memspace   - Memory space of the data to be set to 0 (HOST or DEVICE)
     *
     * @pre   _j_ < _k_ i.e,, _j_ is smaller than the total number of vectors in multivector.
     */
    template <typename ScalarT, typename IdxT>
    int Vector<ScalarT, IdxT>::setToConst(IdxT j, ScalarT C, memory::MemorySpace memspace)
    {
      using namespace memory;
      switch (memspace)
      {
      case HOST:
        if (h_data_ == nullptr)
        {
          out::error() << "Vector::setToConst - host data not allocated\n";
          return 1;
        }
        mem_.setArrayToConstOnHost(&h_data_[n_size_ * j], C, n_size_);
        cpu_updated_[j] = true;
        gpu_updated_[j] = false;
        break;
      case DEVICE:
        if (d_data_ == nullptr)
        {
          out::error() << "Vector::setToConst - device data not allocated\n";
          return 1;
        }
        mem_.setArrayToConstOnDevice(&d_data_[n_size_ * j], C, n_size_);
        cpu_updated_[j] = false;
        gpu_updated_[j] = true;
        break;
      }
      return 0;
    }

    /**
     * @brief Resize vector to `new_n_size`.
     *
     * Use for vectors and multivectors that change size throughout computation.
     *
     * @note Vector needs to have capacity set to maximum expected size.
     * @warning This method is not to be used in vectors who do not own their
     * data.
     *
     * @param[in] new_n_size - New vector length
     *
     * @return 0 if successful, 1 otherwise.
     *
     * @pre `new_n_size` <= `n_capacity_`
     */
    template <typename ScalarT, typename IdxT>
    int Vector<ScalarT, IdxT>::resize(IdxT new_n_size)
    {
      assert(owns_cpu_data_ && owns_gpu_data_
             && "Cannot resize if vector is not owning the data.");

      if (new_n_size > n_capacity_)
      {
        out::error() << "Trying to resize vector to " << new_n_size
                     << " elements but memory allocated only for " << n_capacity_
                     << " elements.\n";
        return 1;
      }
      else
      {
        n_size_ = new_n_size;
        return 0;
      }
    }

    /**
     * @brief copy HOST or DEVICE data of a specified vector in a multivector to _dest_.
     *
     * This function allows to copy data between different memory spaces in one call.
     * For example, you can copy data of vector _i_ from HOST to DEVICE, or from DEVICE to HOST.
     *
     * @param[out] dest      - Pointer to the memory to which data is copied
     * @param[in] i          - Index of a vector in a multivector
     * @param[in] memspaceInSrc   - Memory space (HOST or DEVICE) of the data to be copied
     * @param[in] memspaceOutDst  - Memory space (HOST or DEVICE) to which data is copied
     *
     * @return 0 if successful, 1 otherwise.
     *
     * @pre _i_ < _k_ i.e,, _i_ is smaller than the total number of vectors in multivector.
     * @pre _dest_ is allocated, and the size of _dest_ is at least _n_ (length of a single vector in the multivector).
     * @pre _dest_ is allocated in memspaceOutDst memory space.
     * @post All elements of the vector _i_ are copied to the array _dest_.
     */
    template <typename ScalarT, typename IdxT>
    int Vector<ScalarT, IdxT>::copyToExternal(ScalarT* dest, IdxT i, memory::MemorySpace memspaceSrc, memory::MemorySpace memspaceDst)
    {
      using namespace memory;
      if (i >= k_)
      {
        out::error() << "Trying to copy data for vector " << i << " in multivector but there are only " << k_ << " vectors!\n";
        return 1;
      }
      if (dest == nullptr)
      {
        out::error() << "Trying to copy data for vector " << i << " in multivector but the destination pointer is not allocated!\n";
        return 1;
      }
      ScalarT* data = getData(i, memspaceSrc);
      if (data == nullptr)
      {
        out::error() << "Trying to copy data for vector " << i << " in multivector but the data is not allocated in the source memory space!\n";
        return 1;
      }
      switch (memspaceSrc)
      {
      case HOST:
        switch (memspaceDst)
        {
        case HOST:
          mem_.copyArrayHostToHost(dest, data, n_size_);
          break;
        case DEVICE:
          mem_.copyArrayHostToDevice(dest, data, n_size_);
          break;
        }
        break;
      case DEVICE:
        switch (memspaceDst)
        {
        case HOST:
          mem_.copyArrayDeviceToHost(dest, data, n_size_);
          break;
        case DEVICE:
          mem_.copyArrayDeviceToDevice(dest, data, n_size_);
          break;
        }
        break;
      }
      return 0;
    }

    /**
     * @brief copy HOST or DEVICE data of multivector to _dest_.
     *
     * This function allows to copy data between different memory spaces in one call.
     * For example, you can copy data of multivector from HOST to DEVICE, or from DEVICE to HOST.
     *
     * @param[out] dest      - Pointer to the memory to which data is copied
     * @param[in] memspaceSrc   - Memory space (HOST or DEVICE) of the data to be copied
     * @param[in] memspaceDst  - Memory space (HOST or DEVICE) to which data is copied
     *
     * @return 0 if successful, 1 otherwise.
     *
     * @pre _dest_ is allocated, and the size of _dest_ is at least _n_ * _k_ (total length of all vectors in the multivector).
     * @pre _dest_ is allocated in memspaceOutDst memory space.
     * @post All elements of all vectors in multivector are copied to the array _dest_.
     */
    template <typename ScalarT, typename IdxT>
    int Vector<ScalarT, IdxT>::copyToExternal(ScalarT* dest, memory::MemorySpace memspaceSrc, memory::MemorySpace memspaceDst)
    {
      using namespace memory;
      ScalarT* data = this->getData(memspaceSrc);
      // Check that the source data is not null and up to date
      if (data == nullptr)
      {
        out::error() << "Trying to copy data for multivector but the data is not allocated in the source memory space!\n";
        return 1;
      }
      // Check that the destination memory space is allocated
      if (dest == nullptr)
      {
        out::error() << "Trying to copy data for multivector but the destination pointer is not allocated!\n";
        return 1;
      }
      switch (memspaceSrc)
      {
      case HOST:
        if (!cpu_updated_[0])
        {
          out::error() << "Trying to copy data for multivector but the data is not up to date in the source memory space!\n";
          return 1;
        }
        switch (memspaceDst)
        {
        case HOST:
          mem_.copyArrayHostToHost(dest, data, n_size_ * k_);
          break;
        case DEVICE:
          mem_.copyArrayHostToDevice(dest, data, n_size_ * k_);
          break;
        }
        break;
      case DEVICE:
        if (!gpu_updated_[0])
        {
          out::error() << "Trying to copy data for multivector but the data is not up to date in the source memory space!\n";
          return 1;
        }
        switch (memspaceDst)
        {
        case HOST:
          mem_.copyArrayDeviceToHost(dest, data, n_size_ * k_);
          break;
        case DEVICE:
          mem_.copyArrayDeviceToDevice(dest, data, n_size_ * k_);
          break;
        }
        break;
      }
      return 0;
    }

    //
    // Private methods
    //

    template <typename ScalarT, typename IdxT>
    void Vector<ScalarT, IdxT>::setHostUpdated(bool is_updated)
    {
      std::fill(cpu_updated_, cpu_updated_ + k_, is_updated);
    }

    template <typename ScalarT, typename IdxT>
    void Vector<ScalarT, IdxT>::setDeviceUpdated(bool is_updated)
    {
      std::fill(gpu_updated_, gpu_updated_ + k_, is_updated);
    }

    // template class Vector<double, long int>;
    template class Vector<double, size_t>;
    // template class Vector<double, int>;

  } // namespace LinearAlgebra
} // namespace GridKit
