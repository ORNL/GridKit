gridkit_add_library(
  optimal_power_flow_signal
  INTERFACE_TARGET
  HEADERS Signal.hpp SignalImpl.hpp
  LINK_LIBRARIES INTERFACE GridKit::optimal_power_flow_core)

# Link to interface target for all components
target_link_libraries(
  optimal_power_flow_components
  INTERFACE GridKit::optimal_power_flow_signal)
target_link_libraries(
  optimal_power_flow_components_dependency_tracking
  INTERFACE GridKit::optimal_power_flow_signal)
