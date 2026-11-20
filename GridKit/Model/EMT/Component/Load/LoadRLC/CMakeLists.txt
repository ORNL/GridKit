gridkit_add_library(
  emt_load_rlc
  SOURCES LoadRLCDependencyTracking.cpp
  HEADERS LoadRLC.hpp LoadRLCData.hpp
  LINK_LIBRARIES PUBLIC GridKit::emt_core)

target_link_libraries(
  emt_components
  INTERFACE GridKit::emt_load_rlc)
