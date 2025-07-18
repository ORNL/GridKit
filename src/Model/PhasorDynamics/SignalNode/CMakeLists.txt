
# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#    - Slaven Peles <peless@ornl.gov>
# ]]

gridkit_add_library(phasor_dynamics_signal
  SOURCES
    SignalNode.cpp
  OUTPUT_NAME
    gridkit_phasor_dynamics_signal)

# Link to interface target for all components
target_link_libraries(gridkit_phasor_dynamics_components INTERFACE GRIDKIT::phasor_dynamics_signal)

gridkit_add_library(phasor_dynamics_signal_dependency_tracking
  SOURCES
    SignalNodeDependencyTracking.cpp
  OUTPUT_NAME
    gridkit_phasor_dynamics_signal_dependency_tracking)
