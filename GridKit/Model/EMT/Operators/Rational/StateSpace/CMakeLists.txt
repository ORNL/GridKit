gridkit_add_library(
  emt_state_space
  SOURCES StateSpace.cpp
  HEADERS StateSpace.hpp StateSpaceData.hpp StateSpaceDataJSONParser.hpp
  INCLUDE_DIRECTORIES PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include
  LINK_LIBRARIES PUBLIC GridKit::emt_core GridKit::emt_signal)

gridkit_add_library(
  emt_state_space_dependency_tracking
  SOURCES StateSpaceDependencyTracking.cpp
  INCLUDE_DIRECTORIES PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include
  LINK_LIBRARIES PUBLIC GridKit::emt_core GridKit::emt_signal_dependency_tracking)

target_link_libraries(
  emt_components
  INTERFACE GridKit::emt_state_space)
target_link_libraries(
  emt_components_dependency_tracking
  INTERFACE GridKit::emt_state_space_dependency_tracking)
