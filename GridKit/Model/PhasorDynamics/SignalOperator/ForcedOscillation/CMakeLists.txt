# [[
#  Author(s):
#    - Luke Lowery <lukel@tamu.edu>
# ]]

set(_install_headers ForcedOscillation.hpp ForcedOscillationData.hpp)

if(GRIDKIT_ENABLE_ENZYME)
  gridkit_add_library(
    phasor_dynamics_signal_operator_forced_oscillation
    SOURCES ForcedOscillationEnzyme.cpp
    HEADERS ${_install_headers}
    INCLUDE_DIRECTORIES PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include
    LINK_LIBRARIES
      PUBLIC
      GridKit::phasor_dynamics_core
      PUBLIC
      GridKit::phasor_dynamics_signal
      PRIVATE
      ClangEnzymeFlags
    COMPILE_OPTIONS
      PRIVATE
      -mllvm
      -enzyme-auto-sparsity=1
      -fno-math-errno)
else()
  gridkit_add_library(
    phasor_dynamics_signal_operator_forced_oscillation
    SOURCES ForcedOscillation.cpp
    HEADERS ${_install_headers}
    INCLUDE_DIRECTORIES PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include
    LINK_LIBRARIES
      PUBLIC
      GridKit::phasor_dynamics_core
      PUBLIC
      GridKit::phasor_dynamics_signal)
endif()

gridkit_add_library(
  phasor_dynamics_signal_operator_forced_oscillation_dependency_tracking
  SOURCES ForcedOscillationDependencyTracking.cpp
  INCLUDE_DIRECTORIES PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include
  LINK_LIBRARIES
    PUBLIC
    GridKit::phasor_dynamics_core
    PUBLIC
    GridKit::phasor_dynamics_signal_dependency_tracking)

target_link_libraries(
  phasor_dynamics_components
  INTERFACE GridKit::phasor_dynamics_signal_operator_forced_oscillation)
target_link_libraries(
  phasor_dynamics_components_dependency_tracking
  INTERFACE GridKit::phasor_dynamics_signal_operator_forced_oscillation_dependency_tracking)
