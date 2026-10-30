# [[
#  Author(s):
#    - Luke Lowery <lukel@tamu.edu>
# ]]

set(_install_headers GasCompressor.hpp GasCompressorData.hpp)

if(GRIDKIT_ENABLE_ENZYME)
  set(_gas_compressor_sources GasCompressorEnzyme.cpp)
else()
  set(_gas_compressor_sources GasCompressor.cpp)
endif()

gridkit_add_library(
  phasor_dynamics_gas_compressor
  SOURCES ${_gas_compressor_sources}
  HEADERS ${_install_headers}
  LINK_LIBRARIES PUBLIC GridKit::phasor_dynamics_core)

gridkit_add_library(
  phasor_dynamics_gas_compressor_dependency_tracking
  SOURCES GasCompressorDependencyTracking.cpp
  LINK_LIBRARIES PUBLIC GridKit::phasor_dynamics_core)

# Link to interface target for all components
target_link_libraries(
  phasor_dynamics_components
  INTERFACE GridKit::phasor_dynamics_gas_compressor)
target_link_libraries(
  phasor_dynamics_components_dependency_tracking
  INTERFACE GridKit::phasor_dynamics_gas_compressor_dependency_tracking)
