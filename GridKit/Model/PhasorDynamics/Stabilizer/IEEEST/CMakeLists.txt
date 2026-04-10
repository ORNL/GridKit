# [[
#  Author(s):
#    - Luke Lowery <lukel@tamu.edu>
# ]]

set(_install_headers
  Ieeest.hpp
  IeeestData.hpp)

if(GRIDKIT_ENABLE_ENZYME)
  gridkit_add_library(phasor_dynamics_stabilizer_ieeest
    SOURCES
      IeeestEnzyme.cpp
    HEADERS
      ${_install_headers}
    INCLUDE_DIRECTORIES
      PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include
    LINK_LIBRARIES
      PUBLIC GridKit::phasor_dynamics_core
      PUBLIC GridKit::phasor_dynamics_signal
      PRIVATE ClangEnzymeFlags
    COMPILE_OPTIONS
      PRIVATE -mllvm -enzyme-auto-sparsity=1 -fno-math-errno)
else()
  gridkit_add_library(phasor_dynamics_stabilizer_ieeest
    SOURCES
      Ieeest.cpp
    HEADERS
      ${_install_headers}
    INCLUDE_DIRECTORIES
      PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include
    LINK_LIBRARIES
      PUBLIC GridKit::phasor_dynamics_core
      PUBLIC GridKit::phasor_dynamics_signal)
endif()

gridkit_add_library(phasor_dynamics_stabilizer_ieeest_dependency_tracking
  SOURCES
    IeeestDependencyTracking.cpp
  INCLUDE_DIRECTORIES
    PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include
  LINK_LIBRARIES
    PUBLIC GridKit::phasor_dynamics_core
    PUBLIC GridKit::phasor_dynamics_signal_dependency_tracking)

# Link to interface target for all components
target_link_libraries(phasor_dynamics_components
  INTERFACE GridKit::phasor_dynamics_stabilizer_ieeest)
target_link_libraries(phasor_dynamics_components_dependency_tracking
  INTERFACE GridKit::phasor_dynamics_stabilizer_ieeest_dependency_tracking)
