gridkit_add_library(
  emt_propagation
  SOURCES Propagation.cpp
  HEADERS Propagation.hpp PropagationData.hpp PropagationDataJSONParser.hpp
  INCLUDE_DIRECTORIES PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include
  LINK_LIBRARIES
    PUBLIC
    GridKit::emt_core
    GridKit::emt_delay
    GridKit::emt_vector_fit)

target_link_libraries(
  emt_components
  INTERFACE GridKit::emt_propagation)

gridkit_add_library(
  emt_propagation_dependency_tracking
  SOURCES PropagationDependencyTracking.cpp
  INCLUDE_DIRECTORIES PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include
  LINK_LIBRARIES
    PUBLIC
    GridKit::emt_core
    GridKit::emt_delay_dependency_tracking
    GridKit::emt_vector_fit_dependency_tracking)

target_link_libraries(
  emt_components_dependency_tracking
  INTERFACE GridKit::emt_propagation_dependency_tracking)
