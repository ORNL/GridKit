
if(GRIDKIT_ENABLE_ENZYME)
  gridkit_add_library(phasor_dynamics_bus_fault
    SOURCES
      BusFaultEnzyme.cpp
    LINK_LIBRARIES
      ClangEnzymeFlags
    COMPILE_OPTIONS
      -mllvm -enzyme-auto-sparsity=1 -fno-math-errno)
else()
  gridkit_add_library(phasor_dynamics_bus_fault SOURCES BusFault.cpp)
endif()

gridkit_add_library(phasor_dynamics_bus_fault_dependency_tracking
  SOURCES BusFaultDependencyTracking.cpp)

# Link to interface target for all components
target_link_libraries(gridkit_phasor_dynamics_components
  INTERFACE GRIDKIT::phasor_dynamics_bus_fault)
target_link_libraries(gridkit_phasor_dynamics_components_dependency_tracking
  INTERFACE GRIDKIT::phasor_dynamics_bus_fault_dependency_tracking)
