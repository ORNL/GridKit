# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#    - Slaven Peles <peless@ornl.gov>
# ]]

gridkit_add_library(
  phasor_dynamics_signal
  SOURCES SignalNode.cpp SignalNodeJunction.cpp
  HEADERS SignalNode.hpp SignalNodeData.hpp SignalNodeJunction.hpp
  LINK_LIBRARIES GridKit::phasor_dynamics_core)

gridkit_add_library(
  phasor_dynamics_signal_dependency_tracking
  SOURCES SignalNodeDependencyTracking.cpp SignalNodeJunctionDependencyTracking.cpp
  LINK_LIBRARIES GridKit::phasor_dynamics_core)

# Link to interface target for all components
target_link_libraries(
  phasor_dynamics_components
  INTERFACE GridKit::phasor_dynamics_signal)
target_link_libraries(
  phasor_dynamics_components_dependency_tracking
  INTERFACE GridKit::phasor_dynamics_signal_dependency_tracking)
