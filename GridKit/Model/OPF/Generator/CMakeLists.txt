gridkit_add_library(
  model_opf_generator
  SOURCES GeneratorEnzyme.cpp
  HEADERS Generator.hpp GeneratorData.hpp
  LINK_LIBRARIES
    PUBLIC
    GridKit::model_opf_core
    PRIVATE
    GridKit::optimization_enzyme)

target_link_libraries(
  model_opf_components
  INTERFACE GridKit::model_opf_generator)
