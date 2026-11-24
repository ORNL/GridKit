gridkit_add_library(
  emt_transformer_three_winding
  SOURCES TransformerThreeWindingEnzyme.cpp
  HEADERS TransformerThreeWinding.hpp TransformerThreeWindingData.hpp
  LINK_LIBRARIES
    PRIVATE
    ClangEnzymeFlags
    PUBLIC
    GridKit::emt_core
  COMPILE_OPTIONS
    PRIVATE
    -mllvm
    -enzyme-auto-sparsity=1
    -fno-math-errno)

target_link_libraries(
  emt_components
  INTERFACE GridKit::emt_transformer_three_winding)
