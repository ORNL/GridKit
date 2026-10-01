#

# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#]]

gridkit_add_library(
  solvers_opt
  SOURCES DynamicObjective.cpp DynamicConstraint.cpp
  HEADERS DynamicConstraint.hpp DynamicObjective.hpp OptimizationSolver.hpp
  LINK_LIBRARIES
    PUBLIC
    IPOPT
    PUBLIC
    GridKit::solvers_dyn)
