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

# Link to interface target for all components
target_link_libraries(gridkit_phasor_dynamics_components INTERFACE GRIDKIT::phasor_dynamics_bus_fault)
target_link_libraries(gridkit_phasor_dynamics_components_dependency_tracking INTERFACE GRIDKIT::phasor_dynamics_bus_fault_dependency_tracking)
