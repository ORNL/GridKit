#pragma once
#include <string>

#include <GridKit/MemoryUtilities/MemoryUtils.hpp>

namespace GridKit
{
  namespace LinearAlgebra
  {
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

      ScalarT*       getData(memory::MemorySpace memspace = memory::HOST);
      ScalarT*       getData(IdxT i, memory::MemorySpace memspace = memory::HOST);
      const ScalarT* getData(memory::MemorySpace memspace = memory::HOST) const;
      const ScalarT* getData(IdxT i, memory::MemorySpace memspace = memory::HOST) const;

      IdxT getCapacity() const;
      IdxT getSize() const;
      IdxT getNumVectors() const;

      int setDataUpdated(memory::MemorySpace memspace = memory::HOST);
      int setData(ScalarT* data, memory::MemorySpace memspace = memory::HOST);
      int setData(ScalarT* data, IdxT size, memory::MemorySpace memspace = memory::HOST);
      int allocate(memory::MemorySpace memspace = memory::HOST);
      int setToZero(memory::MemorySpace memspace = memory::HOST);
      int setToZero(IdxT i, memory::MemorySpace memspace = memory::HOST);
      int setToConst(ScalarT C, memory::MemorySpace memspace = memory::HOST);
      int setToConst(IdxT i, ScalarT C, memory::MemorySpace memspace = memory::HOST);
      int syncData(memory::MemorySpace memspaceOut = memory::HOST);
      int resize(IdxT new_n_current);
      int copyToExternal(ScalarT*            dest,
                         IdxT                i,
                         memory::MemorySpace memspaceSrc = memory::HOST,
                         memory::MemorySpace memspaceDst = memory::HOST);
      int copyToExternal(ScalarT*            dest,
                         memory::MemorySpace memspaceSrc = memory::HOST,
                         memory::MemorySpace memspaceDst = memory::HOST);

    private:
      void setHostUpdated(bool is_updated);
      void setDeviceUpdated(bool is_updated);

      IdxT     n_capacity_{0};        ///< vector capacity
      IdxT     k_{0};                 ///< number of vectors in multivector
      IdxT     n_size_{0};            ///< actual size of the vector
      ScalarT* d_data_{nullptr};      ///< DEVICE data array
      ScalarT* h_data_{nullptr};      ///< HOST data array
      bool     gpu_updated_{false};   ///< DEVICE data flags (updated or not)
      bool     cpu_updated_{false};   ///< HOST data flags (updated or not)

      bool owns_gpu_data_{true}; ///< data owneship flag for DEVICE data
      bool owns_cpu_data_{true}; ///< data ownership flag for HOST data

      MemoryManager mem_; ///< Device memory manager object
    };
  } // namespace LinearAlgebra
} // namespace GridKit
