# [[
#  Author(s):
#    - Luke Lowery <lukel@tamu.edu>
#]]


gridkit_add_library(phasor_dynamics_bus
  SOURCES
    BusElectric/BusNetwork.cpp
    BusElectric/BusInfinite.cpp
    BusControl/BusSignal.cpp
  OUTPUT_NAME
    gridkit_phasor_dynamics_bus)