gridkit_example_add_file(${CMAKE_SOURCE_DIR}/cases/PhasorDynamics/WECC240/WECC240Relay.case.json)
gridkit_example_add_file(WECC240Relay.solver.json)

add_test(
  NAME WECC240Relay_line_trip
  COMMAND DynamicSimulation WECC240Relay.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
