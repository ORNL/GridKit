gridkit_add_library(
  emt_conductor
  INTERFACE_TARGET
  HEADERS Conductor.hpp ConductorData.hpp
  LINK_LIBRARIES INTERFACE GridKit::emt_core)
