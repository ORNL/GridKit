#

# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#]]

if(GRIDKIT_ENABLE_SUNDIALS_SPARSE)
  gridkit_add_library(
    solvers_dyn
    SOURCES Ida.cpp
    HEADERS Ida.hpp DynamicSolver.hpp
    LINK_LIBRARIES
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
      GridKit::sparse_matrix
      GridKit::dense_vector)
else()
  gridkit_add_library(
    solvers_dyn
    SOURCES Ida.cpp
    HEADERS Ida.hpp DynamicSolver.hpp
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
      GridKit::sparse_matrix
      GridKit::dense_vector)
endif()

add_subdirectory(Native)
