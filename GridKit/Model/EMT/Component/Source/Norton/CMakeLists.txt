gridkit_add_library(
  emt_norton
  INTERFACE_TARGET
  HEADERS Norton.hpp
  LINK_LIBRARIES INTERFACE GridKit::emt_core GridKit::emt_vector_fit)

gridkit_add_library(
  emt_norton_dependency_tracking
  INTERFACE_TARGET
  LINK_LIBRARIES INTERFACE GridKit::emt_core GridKit::emt_vector_fit_dependency_tracking)

target_link_libraries(
  emt_components
  INTERFACE GridKit::emt_norton)
target_link_libraries(
  emt_components_dependency_tracking
  INTERFACE GridKit::emt_norton_dependency_tracking)
