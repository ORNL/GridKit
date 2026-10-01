gridkit_example_add_file(${CMAKE_SOURCE_DIR}/cases/PhasorDynamics/Toy/TwoBusIeeet1.case.json)
gridkit_example_add_file(TwoBusIeeet1.solver.json)
add_subdirectory(reference)

add_test(
  NAME TwoBusIeeet1_validation
  COMMAND DynamicSimulation TwoBusIeeet1.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
