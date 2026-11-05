gridkit_add_library(
  emt_load_z
  SOURCES LoadZDependencyTracking.cpp
  HEADERS LoadZ.hpp LoadZData.hpp
  LINK_LIBRARIES
    PUBLIC
    GridKit::emt_core
    PUBLIC
    GridKit::emt_vector_fit)

target_link_libraries(
  emt_components
  INTERFACE GridKit::emt_load_z)
