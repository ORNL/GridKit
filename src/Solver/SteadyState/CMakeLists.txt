#

# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#]]

gridkit_add_library(solvers_steady
  SOURCES
    Kinsol.cpp
  LINK_LIBRARIES
    PUBLIC SUNDIALS::kinsol)
