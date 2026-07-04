#pragma once

#include <iostream>
#include <type_traits>

namespace GridKit
{
  namespace memory
  {
    enum MemorySpace
    {
      HOST = 0,
      DEVICE
    };
  } // namespace memory
} // namespace GridKit

namespace GridKit
{
  /**
   * @class MemoryUtils
   *
   * @brief Provides basic memory allocation, free and copy functions.
   *
   * This class provides abstractions for memory management functions for
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
      *v = new T[static_cast<std::size_t>(n)];
      return 0;
    }

    template <typename T>
    int deleteOnHost(T*& v)
    {
      delete[] v;
      v = nullptr;
      return 0;
    }

    /**
     * @brief Copy an array from HOST to HOST.
     *
     * Trivially copyable types (plain scalars, etc.) are copied with a
     * byte-wise memcpy. Types that are not trivially copyable (e.g. a type
     * owning a heap-allocated data structure) would have that structure
     * corrupted by a byte-wise copy, so those are copied element-by-element
     * using the type's own copy assignment instead.
     */
    template <typename I, typename T>
    int copyArrayHostToHost(T* dst, const T* src, I n)
    {
      if constexpr (std::is_trivially_copyable_v<T>)
      {
        std::size_t arraysize = static_cast<std::size_t>(n) * sizeof(T);
        memcpy(dst, src, arraysize);
      }
      else
      {
        for (I i = 0; i < n; ++i)
          dst[i] = src[i];
      }
      return 0;
    }

    /**
     * @brief Set an array on HOST to zero.
     *
     * Trivially copyable types are zeroed with a byte-wise memset. Types
     * that are not trivially copyable would have their internal state
     * corrupted by a byte-wise zeroing, so those are reset to a
     * value-initialized (default-constructed) state element-by-element
     * instead.
     */
    template <typename I, typename T>
    int setZeroArrayOnHost(T* v, I n)
    {
      if constexpr (std::is_trivially_copyable_v<T>)
      {
        std::size_t arraysize = static_cast<std::size_t>(n) * sizeof(T);
        memset(v, 0, arraysize);
      }
      else
      {
        for (I i = 0; i < n; ++i)
          v[i] = T{};
      }
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
  }; // class MemoryUtils

} // namespace GridKit

// Device-side member templates are defined here so they are visible (and
// implicitly instantiated on demand for whatever <I, T> a caller needs) in
// every translation unit that uses MemoryUtils, instead of relying on a
// hand-maintained list of explicit instantiations for one backend.
#include <GridKit/MemoryUtilities/MemoryUtils.tpp>

#ifdef GridKit_ENABLE_GPU

// Check if GPU support is enabled in GridKit and set appropriate device memory manager.
#if defined GridKit_ENABLE_CUDA
#include <GridKit/LinearAlgebra/cuda/CudaMemory.hpp>
using MemoryHandler = GridKit::MemoryUtils<GridKit::memory::Cuda>;
#elif defined GridKit_ENABLE_HIP
#include <GridKit/LinearAlgebra/hip/HipMemory.hpp>
using MemoryHandler = GridKit::MemoryUtils<GridKit::memory::Hip>;
#else
#error Unrecognized device, probably bug in CMake configuration
#endif

#else

// If no GPU support is present, set device memory manager to a dummy object.
#include <GridKit/MemoryUtilities/cpu/CpuMemory.hpp>
using MemoryHandler = GridKit::MemoryUtils<GridKit::memory::Cpu>;

#endif
