gridkit_add_library(
  emt_tower
  INTERFACE_TARGET
  HEADERS Tower.hpp TowerData.hpp
  LINK_LIBRARIES
    INTERFACE
    GridKit::emt_core
    INTERFACE
    GridKit::emt_conductor)
