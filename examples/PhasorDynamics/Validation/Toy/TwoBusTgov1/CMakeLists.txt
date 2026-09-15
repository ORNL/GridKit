gridkit_example_add_file(${CMAKE_SOURCE_DIR}/cases/PhasorDynamics/Toy/TwoBusTgov1.case.json)
gridkit_example_add_file(TwoBusTgov1.solver.json)
add_subdirectory(reference)

add_test(
  NAME TwoBusTgov1_validation
  COMMAND DynamicSimulation TwoBusTgov1.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
