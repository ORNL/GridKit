# [[
#  Author(s):
#    - Philip Fackler <facklerpw@ornl.gov>
# ]]

set(_install_headers BusToSignalAdapter.hpp BusToSignalAdapterData.hpp)

gridkit_add_library(
  phasor_dynamics_bustosignaladapter
  SOURCES BusToSignalAdapter.cpp
  HEADERS ${_install_headers}
  LINK_LIBRARIES GridKit::phasor_dynamics_core GridKit::phasor_dynamics_signal)

gridkit_add_library(
  phasor_dynamics_bustosignaladapter_dependency_tracking
  SOURCES BusToSignalAdapterDependencyTracking.cpp
  LINK_LIBRARIES GridKit::phasor_dynamics_core GridKit::phasor_dynamics_signal_dependency_tracking)

# Link to interface target for all components
target_link_libraries(
  phasor_dynamics_components
  INTERFACE GridKit::phasor_dynamics_bustosignaladapter)
target_link_libraries(
  phasor_dynamics_components_dependency_tracking
  INTERFACE GridKit::phasor_dynamics_bustosignaladapter_dependency_tracking)
