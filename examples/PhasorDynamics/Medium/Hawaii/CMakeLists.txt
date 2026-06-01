gridkit_example_add_file(hawaii.json)
gridkit_example_add_file(hawaii.solver.json)

add_test(
  NAME hawaii_pdsim
  COMMAND DynamicSimulation hawaii.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
