gridkit_add_library(
  emt_ggov1
  SOURCES Ggov1Enzyme.cpp
  HEADERS Ggov1.hpp Ggov1Data.hpp
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
  INTERFACE GridKit::emt_ggov1)
