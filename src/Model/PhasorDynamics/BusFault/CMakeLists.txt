gridkit_add_library(phasor_dynamics_bus_fault
  SOURCES
    BusFault.cpp
  OUTPUT_NAME
    gridkit_phasor_dynamics_bus_fault)

gridkit_add_library(phasor_dynamics_bus_fault_dependency_tracking
  SOURCES
    BusFaultDependencyTracking.cpp
  OUTPUT_NAME
    gridkit_phasor_dynamics_bus_fault_dependency_tracking)
