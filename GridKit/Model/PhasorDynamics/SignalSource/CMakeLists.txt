set(_install_headers ConstantSignalSource.hpp ConstantSignalSourceData.hpp)

gridkit_add_library(
  phasor_dynamics_signalsource
  SOURCES ConstantSignalSource.cpp
  HEADERS ${_install_headers}
  LINK_LIBRARIES GridKit::phasor_dynamics_core GridKit::phasor_dynamics_signal)

gridkit_add_library(
  phasor_dynamics_signalsource_dependency_tracking
  SOURCES ConstantSignalSourceDependencyTracking.cpp
  LINK_LIBRARIES GridKit::phasor_dynamics_core GridKit::phasor_dynamics_signal_dependency_tracking)

# Link to interface target for all components
target_link_libraries(
  phasor_dynamics_components
  INTERFACE GridKit::phasor_dynamics_signalsource)
target_link_libraries(
  phasor_dynamics_components_dependency_tracking
  INTERFACE GridKit::phasor_dynamics_signalsource_dependency_tracking)
