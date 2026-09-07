gridkit_add_library(
  emt_delay
  SOURCES Delay.cpp
  HEADERS Delay.hpp DelayData.hpp DelayDataJSONParser.hpp
  INCLUDE_DIRECTORIES PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include
  LINK_LIBRARIES PUBLIC GridKit::emt_core GridKit::emt_signal)

target_link_libraries(
  emt_components
  INTERFACE GridKit::emt_delay)

gridkit_add_library(
  emt_delay_dependency_tracking
  SOURCES DelayDependencyTracking.cpp
  INCLUDE_DIRECTORIES PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include
  LINK_LIBRARIES PUBLIC GridKit::emt_core GridKit::emt_signal_dependency_tracking)

target_link_libraries(
  emt_components_dependency_tracking
  INTERFACE GridKit::emt_delay_dependency_tracking)
