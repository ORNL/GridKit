/**
 * @file MemoryUtils.cpp
 *
 * Explicitly instantiates the non-template members of
 * MemoryUtils<memory::Cpu>. The device-side member templates
 * (allocateArrayOnDevice, copyArrayDeviceToHost, etc.) are defined in
 * MemoryUtils.tpp, which is included directly by MemoryUtils.hpp, so they
 * are implicitly instantiated on demand for whatever <I, T> combination
 * each caller (CooMatrix, CsrMatrix, Vector, ...) actually needs.
 *
 * @author Slaven Peles <peless@ornl.gov>
 */

#include <iostream>

#include <GridKit/MemoryUtilities/MemoryUtils.hpp>
#include <GridKit/MemoryUtilities/cpu/CpuMemory.hpp>

namespace GridKit
{
  template void MemoryUtils<memory::Cpu>::deviceSynchronize();
  template int  MemoryUtils<memory::Cpu>::getLastDeviceError();
  template int  MemoryUtils<memory::Cpu>::deleteOnDevice(void*);
} // namespace GridKit
