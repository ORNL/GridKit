# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#    - Slaven Peles <peless@ornl.gov>
#]]

add_subdirectory(BusInfinite)
add_subdirectory(BusSignal)

gridkit_add_library(phasor_dynamics_bus
  SOURCES
    Bus.cpp
  OUTPUT_NAME
    gridkit_phasor_dynamics_bus)
