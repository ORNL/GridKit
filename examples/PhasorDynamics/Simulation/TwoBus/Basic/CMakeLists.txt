gridkit_example_add_file(TwoBusBasic.case.json)
gridkit_example_add_file(TwoBusBasic.solver.json)
gridkit_example_add_file(TwoBusBasic.ref.csv)

add_test(
  NAME TwoBusBasic
  COMMAND DynamicSimulation TwoBusBasic.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
