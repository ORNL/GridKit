gridkit_example_add_file(${CMAKE_SOURCE_DIR}/cases/PhasorDynamics/IEEE39/IEEE39.case.json)
gridkit_example_add_file(IEEE39.solver.json)

add_test(
  NAME IEEE39_contingency
  COMMAND ContingencyAnalysis IEEE39.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
