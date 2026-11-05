gridkit_add_library(
  emt_bus
  SOURCES Bus.cpp
  HEADERS Bus.hpp BusData.hpp
  LINK_LIBRARIES PUBLIC GridKit::emt_core)

target_link_libraries(
  emt_components
  INTERFACE GridKit::emt_bus)
