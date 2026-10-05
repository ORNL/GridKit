#

# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#]]

if(GRIDKIT_ENABLE_SUNDIALS_SPARSE)
  gridkit_add_library(
    solvers_dyn
    SOURCES Ida.cpp SplittingStep.cpp
    HEADERS Ida.hpp SplittingStep.hpp DynamicSolver.hpp
    LINK_LIBRARIES
      PUBLIC
      SUNDIALS::nvecserial
      PUBLIC
      SUNDIALS::nvecmanyvector
      PUBLIC
      SUNDIALS::idas
      PUBLIC
      SUNDIALS::kinsol
      PUBLIC
      SUNDIALS::arkode # SUNAdaptController_PI
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
    SOURCES Ida.cpp SplittingStep.cpp
    HEADERS Ida.hpp SplittingStep.hpp DynamicSolver.hpp
    LINK_LIBRARIES
      PUBLIC
      SUNDIALS::nvecserial
      PUBLIC
      SUNDIALS::nvecmanyvector
      PUBLIC
      SUNDIALS::idas
      PUBLIC
      SUNDIALS::kinsol
      PUBLIC
      SUNDIALS::arkode # SUNAdaptController_PI
      PUBLIC
      GridKit::definitions
      PUBLIC
      GridKit::utilities_logger
      PUBLIC
      GridKit::sparse_matrix
      GridKit::dense_vector)
endif()

if(GRIDKIT_ENABLE_OPENMP)
  target_link_libraries(solvers_dyn PUBLIC OpenMP::OpenMP_CXX)
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
