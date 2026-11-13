set(_install_headers Regca.hpp RegcaData.hpp)

gridkit_add_library(
  optimal_power_flow_regca
  SOURCES RegcaEnzyme.cpp
  HEADERS ${_install_headers}
  LINK_LIBRARIES
    PRIVATE
    ClangEnzymeFlags
    PUBLIC
    GridKit::optimal_power_flow_core
  COMPILE_OPTIONS
    PRIVATE
    -mllvm
    -enzyme-auto-sparsity=1
    -fno-math-errno)

gridkit_add_library(
  optimal_power_flow_regca_dependency_tracking
  SOURCES RegcaDependencyTracking.cpp
  LINK_LIBRARIES PUBLIC GridKit::optimal_power_flow_core)

# Link to interface target for all components
target_link_libraries(
  optimal_power_flow_components
  INTERFACE GridKit::optimal_power_flow_regca)
target_link_libraries(
  optimal_power_flow_components_dependency_tracking
  INTERFACE GridKit::optimal_power_flow_regca_dependency_tracking)
