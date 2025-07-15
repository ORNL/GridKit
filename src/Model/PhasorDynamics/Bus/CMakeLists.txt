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

# Link to interface target for all components
target_link_libraries(gridkit_phasor_dynamics_components INTERFACE GRIDKIT::phasor_dynamics_bus)
