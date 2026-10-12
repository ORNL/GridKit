if(GRIDKIT_ENABLE_ENZYME)
  set(_source SignalDelayEnzyme.cpp)
else()
  set(_source SignalDelay.cpp)
endif()

gridkit_add_library(
  phasor_dynamics_signal_delay
  SOURCES ${_source}
  HEADERS SignalDelay.hpp SignalDelayData.hpp SignalDelayControl.hpp
  LINK_LIBRARIES GridKit::phasor_dynamics_core GridKit::phasor_dynamics_signal)

gridkit_add_library(
  phasor_dynamics_signal_delay_dependency_tracking
  SOURCES SignalDelayDependencyTracking.cpp
  LINK_LIBRARIES GridKit::phasor_dynamics_core GridKit::phasor_dynamics_signal_dependency_tracking)

target_link_libraries(
  phasor_dynamics_components
  INTERFACE GridKit::phasor_dynamics_signal_delay)
target_link_libraries(
  phasor_dynamics_components_dependency_tracking
  INTERFACE GridKit::phasor_dynamics_signal_delay_dependency_tracking)
