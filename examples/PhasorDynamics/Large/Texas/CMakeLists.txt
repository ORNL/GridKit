add_test(
  NAME texas_pdsim
  COMMAND DynamicSimulation ${CMAKE_CURRENT_SOURCE_DIR}/texas.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
