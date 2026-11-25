gridkit_add_library(
  emt_signal
  SOURCES SignalDependencyTracking.cpp
  HEADERS Signal.hpp SignalData.hpp
  LINK_LIBRARIES PUBLIC GridKit::emt_core)

target_link_libraries(
  emt_components
  INTERFACE GridKit::emt_signal)
