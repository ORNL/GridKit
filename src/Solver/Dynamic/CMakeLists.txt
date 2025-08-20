#

# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#]]

gridkit_add_library(solvers_dyn
  OPTIONAL
  SOURCES
    Ida.cpp
  LINK_LIBRARIES
    PUBLIC SUNDIALS::nvecserial
    PUBLIC SUNDIALS::idas
    PUBLIC SUNDIALS::sunlinsolklu
  OUTPUT_NAME
    gridkit_solvers_dyn)
