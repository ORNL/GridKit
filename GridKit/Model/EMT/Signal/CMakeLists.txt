# [[
#  Author(s):
#    - Luke Lowery <lukel@tamu.edu>
# ]]

gridkit_add_library(
  emt_signal
  SOURCES Signal.cpp
  HEADERS Signal.hpp SignalData.hpp)

gridkit_add_library(emt_signal_dependency_tracking SOURCES SignalDependencyTracking.cpp)

# Link to interface target for all components
target_link_libraries(
  emt_components
  INTERFACE GridKit::emt_signal)
target_link_libraries(
  emt_components_dependency_tracking
  INTERFACE GridKit::emt_signal_dependency_tracking)
