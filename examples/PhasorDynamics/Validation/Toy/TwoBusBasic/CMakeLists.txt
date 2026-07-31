gridkit_example_add_file(${CMAKE_SOURCE_DIR}/cases/PhasorDynamics/Toy/TwoBusBasic.case.json)
gridkit_example_add_file(TwoBusBasic.solver.json)
add_subdirectory(reference)

add_test(
  NAME TwoBusBasic_validation
  COMMAND DynamicSimulation TwoBusBasic.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

gridkit_example_add_file(${CMAKE_SOURCE_DIR}/cases/PhasorDynamics/Toy/TwoBusBasic.state.json)
