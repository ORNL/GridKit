# [[
#  Author(s):
#    - Luke Lowery <lukel@tamu.edu>
#    - Adam Birchfield <abirchfield@tamu.edu>
# ]]

set(_install_headers
  Ieeet1.hpp
  Ieeet1Data.hpp)

if(GRIDKIT_ENABLE_ENZYME)
  gridkit_add_library(phasor_dynamics_exciter_ieeet1
    SOURCES
    Ieeet1Enzyme.cpp
    HEADERS
      ${_install_headers}
    LINK_LIBRARIES
      PUBLIC GridKit::phasor_dynamics_core
      PUBLIC GridKit::phasor_dynamics_signal
      PRIVATE GridKit::phasor_dynamics_utils
      PRIVATE ClangEnzymeFlags
    COMPILE_OPTIONS
      PRIVATE -mllvm -enzyme-auto-sparsity=1 -fno-math-errno)
else()
  gridkit_add_library(phasor_dynamics_exciter_ieeet1
    SOURCES
      Ieeet1.cpp
    HEADERS
      ${_install_headers}
    LINK_LIBRARIES
      PRIVATE GridKit::phasor_dynamics_utils
      PUBLIC GridKit::phasor_dynamics_core
      PUBLIC GridKit::phasor_dynamics_signal)
endif()

gridkit_add_library(phasor_dynamics_exciter_ieeet1_dependency_tracking
  SOURCES
    Ieeet1DependencyTracking.cpp
  LINK_LIBRARIES
    PRIVATE GridKit::phasor_dynamics_utils
    PUBLIC GridKit::phasor_dynamics_core
    PUBLIC GridKit::phasor_dynamics_signal_dependency_tracking)

# Link to interface target for all components
target_link_libraries(phasor_dynamics_components
  INTERFACE GridKit::phasor_dynamics_exciter_ieeet1)
target_link_libraries(phasor_dynamics_components_dependency_tracking
  INTERFACE GridKit::phasor_dynamics_exciter_ieeet1_dependency_tracking)
