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
      PUBLIC SUNDIALS::sunlinsolklu)
else()
  gridkit_add_library(solvers_dyn
    SOURCES
      Ida.cpp
    LINK_LIBRARIES
      PUBLIC SUNDIALS::nvecserial
      PUBLIC SUNDIALS::idas)
endif()
