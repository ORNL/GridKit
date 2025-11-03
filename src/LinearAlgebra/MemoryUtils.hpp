#pragma once

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
         * @return Always return success!
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
         * @return Always return failure!
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
         * @return Always return failure!
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
         * @return Always return failure!
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
         * @return Always return failure!
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
         * @return Always return failure!
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
         * @return Always return failure!
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
         * @return Always return failure!
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

      }; // struct Cpu
    }; // namespace memory

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

    using MemoryHandler = MemoryUtils<memory::Cpu>;
  } // namespace LinearAlgebra
} // namespace GridKit
