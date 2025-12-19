# [[
#  Author(s):
#    - Luke Lowery <lukel@tamu.edu>
#    - Adam Birchfield <abirchfield@tamu.edu>
# ]]

set(_install_headers
  Ieeet1.hpp
  Ieeet1Data.hpp)

gridkit_add_library(phasor_dynamics_exciter_ieeet1
  SOURCES
    Ieeet1.cpp
  HEADERS
    ${_install_headers}
  INCLUDE_DIRECTORIES
    PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include
  LINK_LIBRARIES
    GridKit::phasor_dynamics_core
    GridKit::phasor_dynamics_signal)

gridkit_add_library(phasor_dynamics_exciter_ieeet1_dependency_tracking
  SOURCES
    Ieeet1DependencyTracking.cpp
  INCLUDE_DIRECTORIES
    PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include
  LINK_LIBRARIES
    GridKit::phasor_dynamics_core
    GridKit::phasor_dynamics_signal_dependency_tracking)

# Link to interface target for all components
target_link_libraries(phasor_dynamics_components
  INTERFACE GridKit::phasor_dynamics_exciter_ieeet1)
target_link_libraries(phasor_dynamics_components_dependency_tracking
  INTERFACE GridKit::phasor_dynamics_exciter_ieeet1_dependency_tracking)
