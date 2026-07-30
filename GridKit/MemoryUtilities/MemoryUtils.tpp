/**
 * @file MemoryUtils.tpp
 *
 * Contains implementation of memory utility functions wrappers.
 * All it does it calls vendor specific functions frm an abstract interface.
 *
 * @author Slaven Peles <peless@ornl.gov>
 */

#pragma once

namespace GridKit
{
  template <typename policy>
  void MemoryUtils<policy>::deviceSynchronize()
  {
    Policy::deviceSynchronize();
  }

  template <typename policy>
  int MemoryUtils<policy>::getLastDeviceError()
  {
    return Policy::getLastDeviceError();
  }

  template <typename policy>
  int MemoryUtils<policy>::deleteOnDevice(void* v)
  {
    return Policy::deleteOnDevice(v);
  }

  template <typename policy>
  template <typename I, typename T>
  int MemoryUtils<policy>::allocateArrayOnDevice(T** v, I n)
  {
    return Policy::template allocateArrayOnDevice<I, T>(v, n);
  }

  template <typename policy>
  template <typename I, typename T>
  int MemoryUtils<policy>::allocateBufferOnDevice(T** v, I n)
  {
    return Policy::template allocateBufferOnDevice<I, T>(v, n);
  }

  template <typename policy>
  template <typename I, typename T>
  int MemoryUtils<policy>::setZeroArrayOnDevice(T* v, I n)
  {
    return Policy::template setZeroArrayOnDevice<I, T>(v, n);
  }

  template <typename policy>
  template <typename I, typename T>
  int MemoryUtils<policy>::setArrayToConstOnDevice(T* v, T c, I n)
  {
    return Policy::template setArrayToConstOnDevice<I, T>(v, c, n);
  }

  template <typename policy>
  template <typename I, typename T>
  int MemoryUtils<policy>::copyArrayDeviceToHost(T* dst, const T* src, I n)
  {
    return Policy::template copyArrayDeviceToHost<I, T>(dst, src, n);
  }

  template <typename policy>
  template <typename I, typename T>
  int MemoryUtils<policy>::copyArrayDeviceToDevice(T* dst, const T* src, I n)
  {
    return Policy::template copyArrayDeviceToDevice<I, T>(dst, src, n);
  }

  template <typename policy>
  template <typename I, typename T>
  int MemoryUtils<policy>::copyArrayHostToDevice(T* dst, const T* src, I n)
  {
    return Policy::template copyArrayHostToDevice<I, T>(dst, src, n);
  }

} // namespace GridKit
