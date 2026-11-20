gridkit_add_library(
  emt_reecb
  SOURCES ReecbEnzyme.cpp
  HEADERS Reecb.hpp ReecbData.hpp
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
  INTERFACE GridKit::emt_reecb)
