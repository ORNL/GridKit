gridkit_add_library(phasor_dynamics_genrou
  SOURCES
    Genrou.cpp
  OUTPUT_NAME
    gridkit_phasor_dynamics_genrou)

# Link to interface target for all components
target_link_libraries(gridkit_phasor_dynamics_components INTERFACE GRIDKIT::phasor_dynamics_genrou)
