gridkit_example_add_file(TwoBusTgov1.case.json)
gridkit_example_add_file(TwoBusTgov1.solver.json)
gridkit_example_add_file(TwoBusTgov1.ref.csv)

add_test(
  NAME TwoBusTgov1Json
  COMMAND DynamicSimulation TwoBusTgov1.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
