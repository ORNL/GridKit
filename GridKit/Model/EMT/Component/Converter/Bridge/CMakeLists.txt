gridkit_add_library(
  emt_bridge
  SOURCES BridgeEnzyme.cpp
  HEADERS Bridge.hpp BridgeData.hpp
  LINK_LIBRARIES
    PRIVATE
    ClangEnzymeFlags
    PUBLIC
    GridKit::emt_core
    PUBLIC
    GridKit::emt_pwm
  COMPILE_OPTIONS
    PRIVATE
    -mllvm
    -enzyme-auto-sparsity=1
    -fno-math-errno)

target_link_libraries(
  emt_components
  INTERFACE GridKit::emt_bridge)
