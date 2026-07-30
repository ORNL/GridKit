# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#    - Slaven Peles <peless@ornl.gov>
#]]

set(_install_headers
    Bus.hpp
    BusData.hpp
    BusFactory.hpp
    BusInfinite.hpp)

if(GRIDKIT_ENABLE_ENZYME)
  gridkit_add_library(
    phasor_dynamics_bus
    SOURCES BusEnzyme.cpp BusInfinite.cpp
    HEADERS ${_install_headers}
    INCLUDE_DIRECTORIES PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include
    LINK_LIBRARIES PUBLIC GridKit::phasor_dynamics_core GridKit::phasor_dynamics_signal)
else()
  gridkit_add_library(
    phasor_dynamics_bus
    SOURCES Bus.cpp BusInfinite.cpp
    HEADERS ${_install_headers}
    INCLUDE_DIRECTORIES PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include
    LINK_LIBRARIES PUBLIC GridKit::phasor_dynamics_core GridKit::phasor_dynamics_signal)
endif()

gridkit_add_library(
  phasor_dynamics_bus_dependency_tracking
  SOURCES BusDependencyTracking.cpp BusInfiniteDependencyTracking.cpp
  INCLUDE_DIRECTORIES PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include
  LINK_LIBRARIES PUBLIC GridKit::phasor_dynamics_core
                 GridKit::phasor_dynamics_signal_dependency_tracking)

# Link to interface target for all components
target_link_libraries(
  phasor_dynamics_components
  INTERFACE GridKit::phasor_dynamics_bus)
target_link_libraries(
  phasor_dynamics_components_dependency_tracking
  INTERFACE GridKit::phasor_dynamics_bus_dependency_tracking)
