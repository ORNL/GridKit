gridkit_example_add_file(ThreeBusConstantSource.case.json)
gridkit_example_add_file(ThreeBusConstantSource.solver.json)

add_test(
  NAME ThreeBusConstantSource
  COMMAND DynamicSimulation ThreeBusConstantSource.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
