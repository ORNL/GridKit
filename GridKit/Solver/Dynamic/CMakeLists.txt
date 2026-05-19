#

# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#]]

set(_install_headers
  DynamicSolver.hpp
  Ida.hpp)

set(_io_install_headers
  IdaDiagnostics.hpp)

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

gridkit_add_library(solvers_dyn_io
  SOURCES
    IdaDiagnostics.cpp
  HEADERS
    ${_io_install_headers}
  LINK_LIBRARIES
    PUBLIC GridKit::solvers_dyn
  INCLUDE_DIRECTORIES
    PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/nlohmann-json/include)
