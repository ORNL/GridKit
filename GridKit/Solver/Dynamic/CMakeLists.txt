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

gridkit_add_library(
  solvers_dyn_io
  SOURCES IdaDiagnostics.cpp
  HEADERS IdaDiagnostics.hpp
  LINK_LIBRARIES PUBLIC GridKit::solvers_dyn
  INCLUDE_DIRECTORIES PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/nlohmann-json/include)

add_subdirectory(Native)
