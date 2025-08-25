#

# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#]]

if (GRIDKIT_ENABLE_SUNDIALS_SPARSE)
  gridkit_add_library(solvers_dyn
    SOURCES
      Ida.cpp
    LINK_LIBRARIES
      PUBLIC SUNDIALS::nvecserial
      PUBLIC SUNDIALS::idas
      PUBLIC SUNDIALS::sunlinsolklu
    OUTPUT_NAME
      gridkit_solvers_dyn)
else()
  gridkit_add_library(solvers_dyn
    SOURCES
      Ida.cpp
    LINK_LIBRARIES
      PUBLIC SUNDIALS::nvecserial
      PUBLIC SUNDIALS::idas
    OUTPUT_NAME
      gridkit_solvers_dyn)
endif()
