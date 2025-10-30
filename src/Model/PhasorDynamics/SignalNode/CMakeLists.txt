
# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#    - Slaven Peles <peless@ornl.gov>
# ]]

gridkit_add_library(phasor_dynamics_signal
  SOURCES
    SignalNode.cpp
  HEADERS
    SignalNode.hpp
    SignalNodeData.hpp)

gridkit_add_library(phasor_dynamics_signal_dependency_tracking
  SOURCES SignalNodeDependencyTracking.cpp)

# Link to interface target for all components
target_link_libraries(phasor_dynamics_components
  INTERFACE GridKit::phasor_dynamics_signal)
target_link_libraries(phasor_dynamics_components_dependency_tracking
  INTERFACE GridKit::phasor_dynamics_signal_dependency_tracking)
