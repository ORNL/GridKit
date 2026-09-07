gridkit_add_library(
  emt_dc_link
  SOURCES DcLink.cpp
  HEADERS DcLink.hpp DcLinkData.hpp
  INCLUDE_DIRECTORIES PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include
  LINK_LIBRARIES
    PUBLIC
    GridKit::emt_core
    PUBLIC
    GridKit::emt_signal)

gridkit_add_library(
  emt_dc_link_dependency_tracking
  SOURCES DcLinkDependencyTracking.cpp
  INCLUDE_DIRECTORIES PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include
  LINK_LIBRARIES
    PUBLIC
    GridKit::emt_core
    PUBLIC
    GridKit::emt_signal_dependency_tracking)

target_link_libraries(
  emt_components
  INTERFACE GridKit::emt_dc_link)
target_link_libraries(
  emt_components_dependency_tracking
  INTERFACE GridKit::emt_dc_link_dependency_tracking)
