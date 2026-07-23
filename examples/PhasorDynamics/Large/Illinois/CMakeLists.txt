add_test(
  NAME illinois_pdsim
  COMMAND DynamicSimulation ${CMAKE_CURRENT_SOURCE_DIR}/illinois.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
