gridkit_example_add_file(${CMAKE_SOURCE_DIR}/cases/PhasorDynamics/Toy/ThreeBusBasic.case.json)
gridkit_example_add_file(ThreeBusBasic.solver.json)

add_test(
  NAME ThreeBusBasic_line_switching
  COMMAND DynamicSimulation ThreeBusBasic.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
