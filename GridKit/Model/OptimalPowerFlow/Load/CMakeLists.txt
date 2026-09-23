set(_install_headers Load.hpp LoadData.hpp)

# Constant kernels have no derivatives, so no Enzyme is needed
gridkit_add_library(
  optimal_power_flow_load
  SOURCES Load.cpp
  HEADERS ${_install_headers}
  LINK_LIBRARIES PUBLIC GridKit::optimal_power_flow_core)

gridkit_add_library(
  optimal_power_flow_load_dependency_tracking
  SOURCES LoadDependencyTracking.cpp
  LINK_LIBRARIES PUBLIC GridKit::optimal_power_flow_core)

# Link to interface target for all components
target_link_libraries(
  optimal_power_flow_components
  INTERFACE GridKit::optimal_power_flow_load)
target_link_libraries(
  optimal_power_flow_components_dependency_tracking
  INTERFACE GridKit::optimal_power_flow_load_dependency_tracking)
