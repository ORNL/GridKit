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

if(GRIDKIT_ENABLE_ENZYME)
  add_subdirectory(Ipopt)
endif()
