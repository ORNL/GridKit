


gridkit_add_library(power_elec_disgen
  SOURCES
    DistributedGenerator.cpp
  LINK_LIBRARIES
    PUBLIC ClangEnzymeFlags
  OUTPUT_NAME
    gridkit_power_elec_disgen)
