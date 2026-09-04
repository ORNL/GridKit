# [[
#  Author(s):
#    - Luke Lowery <lukel@tamu.edu>
# ]]

set(_install_headers Tgov1.hpp Tgov1Data.hpp)

if(GRIDKIT_ENABLE_ENZYME)
  gridkit_add_library(
    emt_tgov1
    SOURCES Tgov1Enzyme.cpp
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
    emt_tgov1
    SOURCES Tgov1.cpp
    HEADERS ${_install_headers}
    INCLUDE_DIRECTORIES PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include
    LINK_LIBRARIES
      PUBLIC
      GridKit::emt_core
      PUBLIC
      GridKit::emt_signal)
endif()

gridkit_add_library(
  emt_tgov1_dependency_tracking
  SOURCES Tgov1DependencyTracking.cpp
  INCLUDE_DIRECTORIES PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include
  LINK_LIBRARIES
    PUBLIC
    GridKit::emt_core
    PUBLIC
    GridKit::emt_signal_dependency_tracking)

# Link to interface target for all components
target_link_libraries(
  emt_components
  INTERFACE GridKit::emt_tgov1)
target_link_libraries(
  emt_components_dependency_tracking
  INTERFACE GridKit::emt_tgov1_dependency_tracking)
