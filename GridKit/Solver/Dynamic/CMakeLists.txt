#

# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#]]

set(_install_headers
  DynamicSolver.hpp
  Ida.hpp)

if (GRIDKIT_ENABLE_SUNDIALS_SPARSE)
  gridkit_add_library(solvers_dyn
    SOURCES
      Ida.cpp
    HEADERS
      ${_install_headers}
    LINK_LIBRARIES
      PUBLIC SUNDIALS::nvecserial
      PUBLIC SUNDIALS::idas
      PUBLIC SUNDIALS::sunlinsolklu
      PUBLIC GridKit::definitions
      PUBLIC GridKit::utilities_logger)
else()
  gridkit_add_library(solvers_dyn
    SOURCES
      Ida.cpp
    HEADERS
      ${_install_headers}
    LINK_LIBRARIES
      PUBLIC SUNDIALS::nvecserial
      PUBLIC SUNDIALS::idas
      PUBLIC GridKit::definitions
      PUBLIC GridKit::utilities_logger)
endif()
