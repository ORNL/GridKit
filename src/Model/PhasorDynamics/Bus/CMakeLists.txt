# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#    - Slaven Peles <peless@ornl.gov>
#]]

if(GRIDKIT_ENABLE_ENZYME)
  gridkit_add_library(phasor_dynamics_bus
    SOURCES
      BusEnzyme.cpp
      BusInfinite.cpp)
else()
  gridkit_add_library(phasor_dynamics_bus
    SOURCES
      Bus.cpp
      BusInfinite.cpp)
endif()

gridkit_add_library(phasor_dynamics_bus_dependency_tracking
  SOURCES
    BusDependencyTracking.cpp
    BusInfiniteDependencyTracking.cpp)

# Link to interface target for all components
target_link_libraries(gridkit_phasor_dynamics_components
  INTERFACE GRIDKIT::phasor_dynamics_bus)
target_link_libraries(gridkit_phasor_dynamics_components_dependency_tracking
  INTERFACE GRIDKIT::phasor_dynamics_bus_dependency_tracking)
