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
  LINK_LIBRARIES
    GRIDKIT::phasor_dynamics_signal)

gridkit_add_library(phasor_dynamics_exciter_ieeet1_dependency_tracking
  SOURCES
    Ieeet1DependencyTracking.cpp
  LINK_LIBRARIES
    GRIDKIT::phasor_dynamics_signal_dependency_tracking)

# Link to interface target for all components
target_link_libraries(gridkit_phasor_dynamics_components
  INTERFACE GRIDKIT::phasor_dynamics_exciter_ieeet1)
target_link_libraries(gridkit_phasor_dynamics_components_dependency_tracking
  INTERFACE GRIDKIT::phasor_dynamics_exciter_ieeet1_dependency_tracking)
