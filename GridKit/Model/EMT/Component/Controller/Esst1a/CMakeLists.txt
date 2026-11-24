gridkit_add_library(
  emt_esst1a
  SOURCES Esst1aEnzyme.cpp
  HEADERS Esst1a.hpp Esst1aData.hpp
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
  INTERFACE GridKit::emt_esst1a)
