#pragma once

#include <GridKit/MemoryUtilities/MemoryUtils.hpp>

#include <resolve/MemoryUtils.hpp>

namespace GridKit
{
  namespace memory
  {
    /**
     * @brief Converts a GridKit \ref MemorySpace to its corresponding Re::Solve MemorySpace
     *
     */
    inline ReSolve::memory::MemorySpace memorySpaceAsResolve(MemorySpace memspace)
    {
      switch (memspace)
      {
      case HOST:
        return ReSolve::memory::HOST;
      case DEVICE:
        return ReSolve::memory::DEVICE;
      default:
        throw "Memory space not supported";
      }
    }
  } // namespace memory
} // namespace GridKit
