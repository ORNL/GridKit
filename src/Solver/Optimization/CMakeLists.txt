#

# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#]]

gridkit_add_library(solvers_opt
  SOURCES
    DynamicObjective.cpp
    DynamicConstraint.cpp
  LINK_LIBRARIES
    PUBLIC Ipopt
    PUBLIC GRIDKIT::solvers_dyn
  OUTPUT_NAME
    gridkit_solvers_opt)
