#[[

@brief Build GridKit backend when there is no GPU support

@author Slaven Peles <peless@ornl.gov>

]]

set(GridKit_CPU_SRC MemoryUtils.cpp)

set(GridKit_CPU_HEADER_INSTALL CpuMemory.hpp)

# First create dummy backend
gridkit_add_library(
  cpu_backend
  SOURCES ${GridKit_CPU_SRC}
  HEADERS ${GridKit_CPU_HEADER_INSTALL})
