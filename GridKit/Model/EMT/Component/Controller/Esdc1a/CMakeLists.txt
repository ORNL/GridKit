gridkit_add_library(
  emt_esdc1a
  SOURCES Esdc1aEnzyme.cpp
  HEADERS Esdc1a.hpp Esdc1aData.hpp
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
  INTERFACE GridKit::emt_esdc1a)
