gridkit_example_add_file(${CMAKE_SOURCE_DIR}/cases/PhasorDynamics/Toy/ThreeBusRelay.case.json)
gridkit_example_add_file(ThreeBusRelay.solver.json)

add_test(
  NAME ThreeBusRelay_line_trip
  COMMAND DynamicSimulation ThreeBusRelay.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
