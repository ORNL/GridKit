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

set(_bus_signal_voltage_out_headers BusSignalVoltageOut/BusSignalVoltageOut.hpp
                                    BusSignalVoltageOut/BusSignalVoltageOutData.hpp)

set(_bus_signal_voltage_in_headers BusSignalVoltageIn/BusSignalVoltageIn.hpp
                                   BusSignalVoltageIn/BusSignalVoltageInData.hpp)

if(GRIDKIT_ENABLE_ENZYME)
  set(_bus_sources
      BusEnzyme.cpp
      BusInfinite.cpp
      BusSignalVoltageOut/BusSignalVoltageOutEnzyme.cpp
      BusSignalVoltageIn/BusSignalVoltageIn.cpp)
else()
  set(_bus_sources
      Bus.cpp
      BusInfinite.cpp
      BusSignalVoltageOut/BusSignalVoltageOut.cpp
      BusSignalVoltageIn/BusSignalVoltageIn.cpp)
endif()

gridkit_add_library(
  phasor_dynamics_bus
  SOURCES ${_bus_sources}
  HEADERS ${_install_headers}
  LINK_LIBRARIES
    PUBLIC
    GridKit::phasor_dynamics_core
    PUBLIC
    GridKit::phasor_dynamics_signal)

gridkit_add_library(
  phasor_dynamics_bus_dependency_tracking
  SOURCES BusDependencyTracking.cpp
          BusInfiniteDependencyTracking.cpp
          BusSignalVoltageOut/BusSignalVoltageOutDependencyTracking.cpp
          BusSignalVoltageIn/BusSignalVoltageInDependencyTracking.cpp
  LINK_LIBRARIES
    PUBLIC
    GridKit::phasor_dynamics_core
    PUBLIC
    GridKit::phasor_dynamics_signal_dependency_tracking)

# Headers in subdirectories keep their relative path on install
install(FILES ${_bus_signal_voltage_out_headers}
        DESTINATION include/GridKit/Model/PhasorDynamics/Bus/BusSignalVoltageOut)
install(FILES ${_bus_signal_voltage_in_headers}
        DESTINATION include/GridKit/Model/PhasorDynamics/Bus/BusSignalVoltageIn)

# Link to interface target for all components
target_link_libraries(
  phasor_dynamics_components
  INTERFACE GridKit::phasor_dynamics_bus)
target_link_libraries(
  phasor_dynamics_components_dependency_tracking
  INTERFACE GridKit::phasor_dynamics_bus_dependency_tracking)
