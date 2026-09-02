# [[
#  Author(s):
#    - Luke Lowery <lukel@tamu.edu>
# ]]

gridkit_add_library(
  emt_signal
  SOURCES SignalNode.cpp
  HEADERS SignalNode.hpp SignalNodeData.hpp)

gridkit_add_library(emt_signal_dependency_tracking SOURCES SignalNodeDependencyTracking.cpp)

# Link to interface target for all components
target_link_libraries(
  emt_components
  INTERFACE GridKit::emt_signal)
target_link_libraries(
  emt_components_dependency_tracking
  INTERFACE GridKit::emt_signal_dependency_tracking)
