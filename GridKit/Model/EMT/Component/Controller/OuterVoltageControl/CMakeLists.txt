gridkit_add_library(
  emt_outer_voltage_control
  SOURCES OuterVoltageControl.cpp
  HEADERS OuterVoltageControl.hpp OuterVoltageControlData.hpp OuterVoltageControlDataJSONParser.hpp
  INCLUDE_DIRECTORIES PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include
  LINK_LIBRARIES
    PUBLIC
    GridKit::emt_core
    PUBLIC
    GridKit::emt_signal)

gridkit_add_library(
  emt_outer_voltage_control_dependency_tracking
  SOURCES OuterVoltageControlDependencyTracking.cpp
  INCLUDE_DIRECTORIES PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include
  LINK_LIBRARIES
    PUBLIC
    GridKit::emt_core
    PUBLIC
    GridKit::emt_signal_dependency_tracking)

target_link_libraries(
  emt_components
  INTERFACE GridKit::emt_outer_voltage_control)
target_link_libraries(
  emt_components_dependency_tracking
  INTERFACE GridKit::emt_outer_voltage_control_dependency_tracking)
