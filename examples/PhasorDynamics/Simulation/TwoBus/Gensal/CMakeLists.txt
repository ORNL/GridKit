gridkit_example_add_file(TwoBusGensal.case.json)
gridkit_example_add_file(TwoBusGensal.solver.json)
gridkit_example_add_file(TwoBusGensal.ref.csv)
add_test(
  NAME TwoBusGensal
  COMMAND DynamicSimulation TwoBusGensal.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
