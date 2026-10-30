# [[
#  Author(s):
#    - Luke Lowery <lukel@tamu.edu>
# ]]

gridkit_add_library(
  phasor_dynamics_gas_bus_infinite
  SOURCES GasBusInfinite.cpp
  HEADERS GasBusInfinite.hpp
  LINK_LIBRARIES PUBLIC GridKit::phasor_dynamics_core)

gridkit_add_library(
  phasor_dynamics_gas_bus_infinite_dependency_tracking
  SOURCES GasBusInfiniteDependencyTracking.cpp
  LINK_LIBRARIES PUBLIC GridKit::phasor_dynamics_core)

# Link to interface target for all components
target_link_libraries(
  phasor_dynamics_components
  INTERFACE GridKit::phasor_dynamics_gas_bus_infinite)
target_link_libraries(
  phasor_dynamics_components_dependency_tracking
  INTERFACE GridKit::phasor_dynamics_gas_bus_infinite_dependency_tracking)
