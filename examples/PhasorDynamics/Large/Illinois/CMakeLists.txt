gridkit_example_add_file(illinois.json)
gridkit_example_add_file(illinois.solver.json)

add_test(
  NAME illinois_pdsim
  COMMAND DynamicSimulation illinois.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
