# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#    - Slaven Peles <peless@ornl.gov>
#]]

gridkit_add_library(phasor_dynamics_bus
  SOURCES
    Bus.cpp
    BusInfinite.cpp
  OUTPUT_NAME
    gridkit_phasor_dynamics_bus)

