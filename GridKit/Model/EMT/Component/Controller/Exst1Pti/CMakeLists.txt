gridkit_add_library(
  emt_exst1_pti
  SOURCES Exst1PtiEnzyme.cpp
  HEADERS Exst1Pti.hpp Exst1PtiData.hpp
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
  INTERFACE GridKit::emt_exst1_pti)
