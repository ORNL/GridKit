gridkit_add_library(phasor_dynamics_genrou
  SOURCES
    Genrou.cpp
  OUTPUT_NAME
    gridkit_phasor_dynamics_genrou)

gridkit_add_library(phasor_dynamics_genrou_dependency_tracking
  SOURCES
    GenrouDependencyTracking.cpp
  OUTPUT_NAME
    gridkit_phasor_dynamics_genrou_dependency_tracking)
