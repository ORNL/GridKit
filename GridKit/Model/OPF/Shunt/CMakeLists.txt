gridkit_add_library(
  model_opf_shunt
  SOURCES ShuntEnzyme.cpp
  HEADERS Shunt.hpp ShuntData.hpp
  LINK_LIBRARIES
    PUBLIC
    GridKit::model_opf_core
    PRIVATE
    GridKit::optimization_enzyme)

target_link_libraries(
  model_opf_components
  INTERFACE GridKit::model_opf_shunt)
