add_test(
  NAME wecc_pdsim
  COMMAND DynamicSimulation ${CMAKE_CURRENT_SOURCE_DIR}/wecc.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
