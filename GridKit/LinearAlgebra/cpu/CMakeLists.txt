#[[

@brief Build GridKit backend when there is no GPU support

@author Slaven Peles <peless@ornl.gov>

]]

set(GridKit_CPU_SRC MemoryUtils.cpp)

set(GridKit_CPU_HEADER_INSTALL CpuMemory.hpp)

# First create dummy backend
gridkit_add_library(
  cpu_backend
  SOURCES MemoryUtils.cpp
  HEADERS CpuMemory.hpp)

# add_library(resolve_backend_cpu SHARED ${GridKit_CPU_SRC})
# target_link_libraries(resolve_backend_cpu PRIVATE resolve_logger)

# install include headers
# install(FILES ${GridKit_CPU_HEADER_INSTALL} DESTINATION include/resolve/cpu)
