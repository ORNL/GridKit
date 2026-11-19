gridkit_add_library(
  emt_tgov1
  SOURCES Tgov1Enzyme.cpp
  HEADERS Tgov1.hpp Tgov1Data.hpp
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
  INTERFACE GridKit::emt_tgov1)
