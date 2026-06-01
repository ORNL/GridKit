add_library(power_electronics_circuit_node INTERFACE)
target_include_directories(
  power_electronics_circuit_node
  INTERFACE $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}> $<INSTALL_INTERFACE:include>)

add_library(GridKit::power_electronics_circuit_node ALIAS power_electronics_circuit_node)

add_subdirectory(Bus)
add_subdirectory(Capacitor)
add_subdirectory(Resistor)
add_subdirectory(VoltageSource)
add_subdirectory(Inductor)
add_subdirectory(LinearTransformer)
add_subdirectory(InductionMotor)
add_subdirectory(SynchronousMachine)
add_subdirectory(DistributedGenerator)
add_subdirectory(TransmissionLine)
add_subdirectory(MicrogridLoad)
add_subdirectory(MicrogridLine)
add_subdirectory(MicrogridBusDQ)

install(
  FILES CircuitComponent.hpp
        CircuitNode.hpp
        CircuitGraph.hpp
        SystemModelPowerElectronics.hpp
        NodeBase.hpp
  DESTINATION include/GridKit/Model/PowerElectronics)
