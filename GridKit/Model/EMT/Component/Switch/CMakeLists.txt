gridkit_add_library(
  emt_switch
  SOURCES SwitchDependencyTracking.cpp
  HEADERS Switch.hpp SwitchData.hpp
  LINK_LIBRARIES PUBLIC GridKit::emt_core)

target_link_libraries(
  emt_components
  INTERFACE GridKit::emt_switch)
