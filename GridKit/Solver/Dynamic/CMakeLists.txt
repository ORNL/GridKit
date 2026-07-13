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
      GridKit::sparse_matrix)
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
      GridKit::sparse_matrix)
endif()

if(GRIDKIT_ENABLE_RESOLVE)
  gridkit_add_library(
    rosenbrock
    SOURCES Rosenbrock.cpp RosenbrockTableaus.cpp
    HEADERS Rosenbrock.hpp
    LINK_LIBRARIES
      PRIVATE
      ReSolve::ReSolve
      PUBLIC
      GridKit::dense_vector
      PUBLIC
      GridKit::definitions
      PUBLIC
      GridKit::utilities_logger
      PUBLIC
      GridKit::sparse_matrix)

  gridkit_add_library(
    inf_norm
    SOURCES InfNorm.cpp
    HEADERS InfNorm.hpp
    LINK_LIBRARIES
      PRIVATE
      ReSolve::ReSolve
      PUBLIC
      GridKit::dense_vector)
endif()

gridkit_add_library(
  adaptive_step
  SOURCES AdaptiveStep.cpp
  HEADERS AdaptiveStep.hpp)

gridkit_add_library(
  fixed_step
  SOURCES FixedStep.cpp
  HEADERS FixedStep.hpp)
