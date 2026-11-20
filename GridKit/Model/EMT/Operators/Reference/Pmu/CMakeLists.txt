gridkit_add_library(
  emt_pmu
  SOURCES PmuEnzyme.cpp
  HEADERS Pmu.hpp PmuData.hpp
  LINK_LIBRARIES
    PRIVATE
    ClangEnzymeFlags
    PUBLIC
    GridKit::emt_core
    PUBLIC
    GridKit::emt_pll
  COMPILE_OPTIONS
    PRIVATE
    -mllvm
    -enzyme-auto-sparsity=1
    -fno-math-errno)

target_link_libraries(
  emt_components
  INTERFACE GridKit::emt_pmu)
