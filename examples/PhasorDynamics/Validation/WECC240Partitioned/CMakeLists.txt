gridkit_example_add_file(${CMAKE_SOURCE_DIR}/cases/PhasorDynamics/WECC240/WECC240.case.json)
gridkit_example_add_file(${CMAKE_SOURCE_DIR}/cases/PhasorDynamics/WECC240/WECC240.partition.json)
gridkit_example_add_file(WECC240Partitioned.solver.json)

add_test(
  NAME WECC240Partitioned_validation
  COMMAND DynamicSimulation WECC240Partitioned.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
