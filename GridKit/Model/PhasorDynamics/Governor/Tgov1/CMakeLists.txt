# [[
#  Author(s):
#    - Luke Lowery <lukel@tamu.edu>
#    - Adam Birchfield <abirchfield@tamu.edu>
# ]]

set(_install_headers Tgov1.hpp Tgov1Data.hpp)

if(GRIDKIT_ENABLE_ENZYME)
  gridkit_add_library(
    phasor_dynamics_governor_tgov1
    SOURCES Tgov1Enzyme.cpp
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
    phasor_dynamics_governor_tgov1
    SOURCES Tgov1.cpp
    HEADERS ${_install_headers}
    LINK_LIBRARIES GridKit::phasor_dynamics_core GridKit::phasor_dynamics_signal)
endif()

gridkit_add_library(
  phasor_dynamics_governor_tgov1_dependency_tracking
  SOURCES Tgov1DependencyTracking.cpp
  LINK_LIBRARIES
    PUBLIC
    GridKit::phasor_dynamics_core
    PUBLIC
    GridKit::phasor_dynamics_signal_dependency_tracking)

# Link to interface target for all components
target_link_libraries(
  phasor_dynamics_components
  INTERFACE GridKit::phasor_dynamics_governor_tgov1)
target_link_libraries(
  phasor_dynamics_components_dependency_tracking
  INTERFACE GridKit::phasor_dynamics_governor_tgov1_dependency_tracking)
