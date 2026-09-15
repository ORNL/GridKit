gridkit_example_add_file(
  ${CMAKE_SOURCE_DIR}/cases/PhasorDynamics/Toy/ThreeBusConstantSource.case.json)
gridkit_example_add_file(ThreeBusConstantSource.solver.json)

add_test(
  NAME ThreeBusConstantSource
  COMMAND DynamicSimulation ThreeBusConstantSource.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
