#

# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#]]

set(_install_headers DynamicSolver.hpp Ida.hpp SundialsException.hpp)

if(GRIDKIT_ENABLE_SUNDIALS_SPARSE)
  set(_sources Ida.cpp)
  set(_link_libraries
      PUBLIC
      SUNDIALS::nvecserial
      PUBLIC
      SUNDIALS::idas
      PUBLIC
      SUNDIALS::sunlinsolklu
      PUBLIC
      GridKit::definitions
      PUBLIC
      GridKit::utilities_logger
      PUBLIC
      GridKit::sparse_matrix)

  if(TARGET SUNDIALS::arkode)
    list(
      APPEND
      _sources
      Arkode.cpp
      ArkStep.cpp
      ErkStep.cpp)
    list(
      APPEND
      _install_headers
      Arkode.hpp
      ArkStep.hpp
      ErkStep.hpp)
    list(
      APPEND
      _link_libraries
      PUBLIC
      SUNDIALS::arkode)
  endif()

  gridkit_add_library(
    solvers_dyn
    SOURCES ${_sources}
    HEADERS ${_install_headers}
    LINK_LIBRARIES ${_link_libraries})
else()
  gridkit_add_library(
    solvers_dyn
    SOURCES Ida.cpp
    HEADERS ${_install_headers}
    LINK_LIBRARIES
      PUBLIC
      SUNDIALS::nvecserial
      PUBLIC
      SUNDIALS::idas
      PUBLIC
      GridKit::definitions
      PUBLIC
      GridKit::utilities_logger
      PUBLIC
      GridKit::sparse_matrix)
endif()
