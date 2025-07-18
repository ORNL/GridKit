gridkit_add_library(phasor_dynamics_genrou
  SOURCES
    Genrou.cpp
  LINK_LIBRARIES
    PUBLIC GRIDKIT::phasor_dynamics_signal
  OUTPUT_NAME
    gridkit_phasor_dynamics_genrou)

gridkit_add_library(phasor_dynamics_genrou_dependency_tracking
  SOURCES
    GenrouDependencyTracking.cpp
  LINK_LIBRARIES
    PUBLIC GRIDKIT::phasor_dynamics_signal_dependency_tracking
  OUTPUT_NAME
    gridkit_phasor_dynamics_genrou_dependency_tracking)

# Link to interface target for all components
target_link_libraries(gridkit_phasor_dynamics_components INTERFACE GRIDKIT::phasor_dynamics_genrou)
