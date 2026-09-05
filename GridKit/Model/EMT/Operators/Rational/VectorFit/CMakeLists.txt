# The real rational realization supplies an exact Jacobian in every build.
gridkit_add_library(
  emt_vector_fit
  SOURCES VectorFit.cpp
  HEADERS VectorFit.hpp VectorFitData.hpp VectorFitDataJSONParser.hpp
  INCLUDE_DIRECTORIES PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include
  LINK_LIBRARIES PUBLIC GridKit::emt_core GridKit::emt_signal)

gridkit_add_library(
  emt_vector_fit_dependency_tracking
  SOURCES VectorFitDependencyTracking.cpp
  INCLUDE_DIRECTORIES PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include
  LINK_LIBRARIES PUBLIC GridKit::emt_core GridKit::emt_signal_dependency_tracking)

target_link_libraries(
  emt_components
  INTERFACE GridKit::emt_vector_fit)
target_link_libraries(
  emt_components_dependency_tracking
  INTERFACE GridKit::emt_vector_fit_dependency_tracking)
