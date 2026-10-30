# [[
#  Author(s):
#    - Luke Lowery <lukel@tamu.edu>
# ]]

set(_install_headers GasBus.hpp GasBusData.hpp GasBusDataJSONParser.hpp)

if(GRIDKIT_ENABLE_ENZYME)
  set(_gas_bus_sources GasBusEnzyme.cpp)
else()
  set(_gas_bus_sources GasBus.cpp)
endif()

gridkit_add_library(
  phasor_dynamics_gas_bus
  SOURCES ${_gas_bus_sources}
  HEADERS ${_install_headers}
  LINK_LIBRARIES PUBLIC GridKit::phasor_dynamics_core)

gridkit_add_library(
  phasor_dynamics_gas_bus_dependency_tracking
  SOURCES GasBusDependencyTracking.cpp
  LINK_LIBRARIES PUBLIC GridKit::phasor_dynamics_core)

# Link to interface target for all components
target_link_libraries(
  phasor_dynamics_components
  INTERFACE GridKit::phasor_dynamics_gas_bus)
target_link_libraries(
  phasor_dynamics_components_dependency_tracking
  INTERFACE GridKit::phasor_dynamics_gas_bus_dependency_tracking)
