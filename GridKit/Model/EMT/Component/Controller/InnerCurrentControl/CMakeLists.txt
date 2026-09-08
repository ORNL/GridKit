gridkit_add_library(
  emt_inner_current_control
  SOURCES InnerCurrentControl.cpp
  HEADERS InnerCurrentControl.hpp InnerCurrentControlData.hpp InnerCurrentControlDataJSONParser.hpp
  INCLUDE_DIRECTORIES PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include
  LINK_LIBRARIES
    PUBLIC
    GridKit::emt_core
    PUBLIC
    GridKit::emt_signal)

gridkit_add_library(
  emt_inner_current_control_dependency_tracking
  SOURCES InnerCurrentControlDependencyTracking.cpp
  INCLUDE_DIRECTORIES PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include
  LINK_LIBRARIES
    PUBLIC
    GridKit::emt_core
    PUBLIC
    GridKit::emt_signal_dependency_tracking)

target_link_libraries(
  emt_components
  INTERFACE GridKit::emt_inner_current_control)
target_link_libraries(
  emt_components_dependency_tracking
  INTERFACE GridKit::emt_inner_current_control_dependency_tracking)
