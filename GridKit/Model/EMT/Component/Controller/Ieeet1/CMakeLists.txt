gridkit_add_library(
  emt_ieeet1
  SOURCES Ieeet1Enzyme.cpp
  HEADERS Ieeet1.hpp Ieeet1Data.hpp
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
  INTERFACE GridKit::emt_ieeet1)
