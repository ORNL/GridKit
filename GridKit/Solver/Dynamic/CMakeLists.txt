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

if(GRIDKIT_ENABLE_RESOLVE)
  gridkit_add_library(
    rosenbrock
    SOURCES Rosenbrock.cpp RosenbrockTableaus.cpp
    HEADERS Rosenbrock.hpp
            Native/ErrorNorm.hpp
            Native/StepControl.hpp
            Native/StepController.hpp
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
endif()

if(GRIDKIT_ENABLE_RESOLVE)
  gridkit_add_library(
    MultiStageCosim
    SOURCES MultiStageCosimulation.cpp
    HEADERS MultiStageCosimulation.hpp
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
      GridKit::sparse_matrix
      PUBLIC
      GridKit::rosenbrock
      PUBLIC
      GridKit::adaptive_step
      PUBLIC
      GridKit::resolve_system_solver
      PUBLIC
      GridKit::inf_norm
      PUBLIC
      solvers_dyn
      PUBLIC
      GridKit::fixed_step)
endif()

if(GRIDKIT_ENABLE_OPENMP)
  target_link_libraries(MultiStageCosim PRIVATE OpenMP::OpenMP_CXX)
endif()

add_subdirectory(Native)
