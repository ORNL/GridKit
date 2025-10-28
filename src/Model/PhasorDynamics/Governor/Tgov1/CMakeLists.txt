# [[
#  Author(s):
#    - Luke Lowery <lukel@tamu.edu>
#    - Adam Birchfield <abirchfield@tamu.edu>
# ]]

if(GRIDKIT_ENABLE_ENZYME)
  gridkit_add_library(phasor_dynamics_governortgov1
    SOURCES
      Tgov1Enzyme.cpp
    LINK_LIBRARIES
      GRIDKIT::phasor_dynamics_signal
      ClangEnzymeFlags
    COMPILE_OPTIONS
      -mllvm -enzyme-auto-sparsity=1 -fno-math-errno)
else()
  gridkit_add_library(phasor_dynamics_governortgov1
    SOURCES
      Tgov1.cpp
    LINK_LIBRARIES
      GRIDKIT::phasor_dynamics_signal)
endif()

gridkit_add_library(phasor_dynamics_governortgov1_dependency_tracking
  SOURCES
    Tgov1DependencyTracking.cpp
  LINK_LIBRARIES
    PUBLIC GRIDKIT::phasor_dynamics_signal_dependency_tracking)

# Link to interface target for all components
target_link_libraries(gridkit_phasor_dynamics_components
  INTERFACE GRIDKIT::phasor_dynamics_governortgov1)
target_link_libraries(gridkit_phasor_dynamics_components_dependency_tracking
  INTERFACE GRIDKIT::phasor_dynamics_governortgov1_dependency_tracking)
