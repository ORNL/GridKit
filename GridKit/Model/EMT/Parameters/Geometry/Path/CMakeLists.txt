gridkit_add_library(
  emt_path
  INTERFACE_TARGET
  HEADERS Path.hpp PathData.hpp
  LINK_LIBRARIES
    INTERFACE
    GridKit::emt_core
    INTERFACE
    GridKit::emt_tower)
