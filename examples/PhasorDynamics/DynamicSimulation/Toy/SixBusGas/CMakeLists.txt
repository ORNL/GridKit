gridkit_example_add_file(${CMAKE_SOURCE_DIR}/cases/PhasorDynamics/Toy/SixBusGas.case.json)
gridkit_example_add_file(SixBusGas.solver.json)

add_test(
  NAME SixBusGas_supply_pressure_drop
  COMMAND DynamicSimulation SixBusGas.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
