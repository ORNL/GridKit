gridkit_add_library(phasor_dynamics_genrou
  SOURCES
    Genrou.cpp
  LINK_LIBRARIES
    GRIDKIT::phasor_dynamics_governortgov1
  OUTPUT_NAME
    gridkit_phasor_dynamics_genrou)
