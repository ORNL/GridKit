gridkit_add_library(
  emt_geometric_inductance
  SOURCES GeometricInductance.cpp
  HEADERS GeometricInductance.hpp
  LINK_LIBRARIES
    PUBLIC
    GridKit::emt_core
    GridKit::emt_tower
    GridKit::emt_conductor)
