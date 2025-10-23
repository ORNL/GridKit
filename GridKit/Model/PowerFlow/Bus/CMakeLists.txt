#

# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#]]

gridkit_add_library(bus
  SOURCES
    BusPQ.cpp
    BusPV.cpp
    BusSlack.cpp
  OUTPUT_NAME
    gridkit_bus)
