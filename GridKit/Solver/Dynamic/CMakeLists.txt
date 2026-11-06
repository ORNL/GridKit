#

# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#]]

set(_dynamic_solver_sources Ida.cpp)
set(_dynamic_solver_headers Ida.hpp DynamicSolver.hpp)
set(_dynamic_solver_libraries
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

if(GRIDKIT_ENABLE_SUNDIALS_SPARSE)
  list(
    APPEND
    _dynamic_solver_libraries
    PUBLIC
    SUNDIALS::sunlinsolklu)
endif()

gridkit_add_library(
  solvers_dyn
  SOURCES ${_dynamic_solver_sources}
  HEADERS ${_dynamic_solver_headers}
  LINK_LIBRARIES ${_dynamic_solver_libraries})

# PDAEStep creates IDA memory, whereas Ida.cpp requires IDAS for adjoints.
# IDA and IDAS export the same IDA* symbols but have incompatible memory layouts.
# Keep their dependencies separate so partition calls bind to the IDA backend.
if(GRIDKIT_ENABLE_SUNDIALS_PARTITIONED_DAE)
  gridkit_add_library(
    solvers_partitioned
    SOURCES Partitioned.cpp
    HEADERS Partitioned.hpp
    LINK_LIBRARIES
      PUBLIC
      SUNDIALS::ida
      PUBLIC
      SUNDIALS::arkode
      PUBLIC
      SUNDIALS::nvecserial
      PUBLIC
      SUNDIALS::nvecmanyvector
      PUBLIC
      SUNDIALS::sunlinsolspgmr
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

add_subdirectory(Native)
