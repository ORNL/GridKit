gridkit_example_add_file(
  ${CMAKE_SOURCE_DIR}/cases/PhasorDynamics/Toy/ThreeBusOutOfService.case.json)
gridkit_example_add_file(
  ${CMAKE_SOURCE_DIR}/cases/PhasorDynamics/Toy/ThreeBusOutOfService.state.json)
gridkit_example_add_file(ThreeBusOutOfService.solver.json)

add_test(
  NAME ThreeBusOutOfService_validation
  COMMAND DynamicSimulation ThreeBusOutOfService.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
