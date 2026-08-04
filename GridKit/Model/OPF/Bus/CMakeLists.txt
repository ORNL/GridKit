gridkit_add_library(
  model_opf_bus
  SOURCES Bus.cpp
  HEADERS Bus.hpp BusData.hpp
  LINK_LIBRARIES PUBLIC GridKit::model_opf_core)

target_link_libraries(
  model_opf_components
  INTERFACE GridKit::model_opf_bus)
