gridkit_example_add_file(TwoBusIeeet1.case.json)
gridkit_example_add_file(TwoBusIeeet1.solver.json)
gridkit_example_add_file(TwoBusIeeet1.ref.csv)

add_test(
  NAME TwoBusIeeet1
  COMMAND DynamicSimulation TwoBusIeeet1.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
