# [[
#  Author(s):
#    - Luke Lowery <lukel@tamu.edu>
# ]]

set(_install_headers Repca.hpp RepcaData.hpp)

if(GRIDKIT_ENABLE_ENZYME)
  gridkit_add_library(
    phasor_dynamics_controller_repca
    SOURCES RepcaEnzyme.cpp
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
    phasor_dynamics_controller_repca
    SOURCES Repca.cpp
    HEADERS ${_install_headers}
    LINK_LIBRARIES
      PUBLIC
      GridKit::phasor_dynamics_core
      PUBLIC
      GridKit::phasor_dynamics_signal)
endif()

gridkit_add_library(
  phasor_dynamics_controller_repca_dependency_tracking
  SOURCES RepcaDependencyTracking.cpp
  LINK_LIBRARIES
    PUBLIC
    GridKit::phasor_dynamics_core
    PUBLIC
    GridKit::phasor_dynamics_signal_dependency_tracking)

target_link_libraries(
  phasor_dynamics_components
  INTERFACE GridKit::phasor_dynamics_controller_repca)
target_link_libraries(
  phasor_dynamics_components_dependency_tracking
  INTERFACE GridKit::phasor_dynamics_controller_repca_dependency_tracking)
