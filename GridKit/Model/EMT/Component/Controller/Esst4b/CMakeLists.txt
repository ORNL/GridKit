gridkit_add_library(
  emt_esst4b
  SOURCES Esst4bEnzyme.cpp
  HEADERS Esst4b.hpp Esst4bData.hpp
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
  INTERFACE GridKit::emt_esst4b)
