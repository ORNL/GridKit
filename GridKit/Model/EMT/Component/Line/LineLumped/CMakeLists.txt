gridkit_add_library(
  emt_line_lumped
  SOURCES LineLumpedDependencyTracking.cpp
  HEADERS LineLumped.hpp LineLumpedData.hpp
  LINK_LIBRARIES
    PUBLIC
    GridKit::emt_core
    PUBLIC
    GridKit::emt_vector_fit)

target_link_libraries(
  emt_components
  INTERFACE GridKit::emt_line_lumped)
