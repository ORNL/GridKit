gridkit_example_add_file(ThreeBusBasic.case.json)
gridkit_example_add_file(ThreeBusBasic.ref.csv)
gridkit_example_add_file(ThreeBusBasic.solver.json)

add_test(
  NAME ThreeBusBasic
  COMMAND DynamicSimulation ThreeBusBasic.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
