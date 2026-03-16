#

# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#]]

set(_install_headers
  DynamicSolver.hpp
  Ida.hpp
  Rosenbrock.hpp)

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
      PUBLIC GridKit::utilities_logger
      PUBLIC GridKit::sparse_matrix)
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
      PUBLIC GridKit::utilities_logger
      PUBLIC GridKit::sparse_matrix)
endif()

if (GRIDKIT_ENABLE_RESOLVE)
  gridkit_add_library(rosenbrock
    SOURCES
      Rosenbrock.cpp
    HEADERS
      ${_install_headers}
    LINK_LIBRARIES
      PUBLIC ReSolve::ReSolve
      PUBLIC GridKit::definitions
      PUBLIC GridKit::utilities_logger
      PUBLIC GridKit::sparse_matrix
    )
endif()