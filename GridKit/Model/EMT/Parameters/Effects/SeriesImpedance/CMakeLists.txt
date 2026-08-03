gridkit_add_library(
  emt_series_impedance
  SOURCES SeriesImpedance.cpp
  HEADERS SeriesImpedance.hpp
  LINK_LIBRARIES
    PUBLIC
    GridKit::emt_core
    GridKit::emt_skin_effect
    GridKit::emt_geometric_inductance
    GridKit::emt_carson)
