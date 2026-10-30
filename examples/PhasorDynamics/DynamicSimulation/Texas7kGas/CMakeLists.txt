gridkit_example_add_file(${CMAKE_SOURCE_DIR}/cases/PhasorDynamics/Texas7kGas/Texas7kGas.case.json)
gridkit_example_add_file(Texas7kGas.solver.json)

add_test(
  NAME Texas7kGas_rupture
  COMMAND DynamicSimulation Texas7kGas.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
