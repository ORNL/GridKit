#

# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#]]

if(TARGET SUNDIALS::idas)
  gridkit_add_library(
    solvers_opt
    SOURCES DynamicObjective.cpp DynamicConstraint.cpp
    HEADERS DynamicConstraint.hpp DynamicObjective.hpp OptimizationSolver.hpp
    LINK_LIBRARIES
      PUBLIC
      IPOPT
      PUBLIC
      GridKit::solvers_dyn)
endif()

gridkit_add_library(
  solvers_opt_problem
  SOURCES OptimizationProblem.cpp
  HEADERS OptimizationProblem.hpp
  LINK_LIBRARIES
    PUBLIC
    IPOPT
    PUBLIC
    GridKit::dense_vector
    PUBLIC
    GridKit::sparse_matrix)
