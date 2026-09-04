# [[
#  Author(s):
#    - Luke Lowery <lukel@tamu.edu>
# ]]

set(_install_headers Bus.hpp BusData.hpp)

if(GRIDKIT_ENABLE_ENZYME)
  gridkit_add_library(
    emt_bus
    SOURCES BusEnzyme.cpp
    HEADERS ${_install_headers}
    INCLUDE_DIRECTORIES PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include
    LINK_LIBRARIES
      PUBLIC
      GridKit::emt_core
      PUBLIC
      GridKit::emt_signal
      PRIVATE
      ClangEnzymeFlags
    COMPILE_OPTIONS
      PRIVATE
      -mllvm
      -enzyme-auto-sparsity=1
      -fno-math-errno
      -fno-vectorize)
else()
  gridkit_add_library(
    emt_bus
    SOURCES Bus.cpp
    HEADERS ${_install_headers}
    INCLUDE_DIRECTORIES PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include
    LINK_LIBRARIES
      PUBLIC
      GridKit::emt_core
      PUBLIC
      GridKit::emt_signal)
endif()

gridkit_add_library(
  emt_bus_dependency_tracking
  SOURCES BusDependencyTracking.cpp
  INCLUDE_DIRECTORIES PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include
  LINK_LIBRARIES
    PUBLIC
    GridKit::emt_core
    PUBLIC
    GridKit::emt_signal_dependency_tracking)

# Link to interface target for all components
target_link_libraries(
  emt_components
  INTERFACE GridKit::emt_bus)
target_link_libraries(
  emt_components_dependency_tracking
  INTERFACE GridKit::emt_bus_dependency_tracking)
