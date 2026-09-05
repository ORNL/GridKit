gridkit_example_add_file(${CMAKE_SOURCE_DIR}/cases/PhasorDynamics/Hawaii/Hawaii.case.json)
gridkit_example_add_file(Hawaii.solver.json)
add_subdirectory(reference)

add_test(
  NAME Hawaii_validation
  COMMAND DynamicSimulation Hawaii.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
