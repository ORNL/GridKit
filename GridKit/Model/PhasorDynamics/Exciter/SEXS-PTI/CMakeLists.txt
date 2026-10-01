# [[
#  Author(s):
#    - Luke Lowery <lukel@tamu.edu>
# ]]

set(_install_headers SexsPti.hpp SexsPtiData.hpp)

if(GRIDKIT_ENABLE_ENZYME)
  gridkit_add_library(
    phasor_dynamics_exciter_sexspti
    SOURCES SexsPtiEnzyme.cpp
    HEADERS ${_install_headers}
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
    phasor_dynamics_exciter_sexspti
    SOURCES SexsPti.cpp
    HEADERS ${_install_headers}
    LINK_LIBRARIES GridKit::phasor_dynamics_core GridKit::phasor_dynamics_signal)
endif()

gridkit_add_library(
  phasor_dynamics_exciter_sexspti_dependency_tracking
  SOURCES SexsPtiDependencyTracking.cpp
  LINK_LIBRARIES GridKit::phasor_dynamics_core GridKit::phasor_dynamics_signal_dependency_tracking)

target_link_libraries(
  phasor_dynamics_components
  INTERFACE GridKit::phasor_dynamics_exciter_sexspti)
target_link_libraries(
  phasor_dynamics_components_dependency_tracking
  INTERFACE GridKit::phasor_dynamics_exciter_sexspti_dependency_tracking)
