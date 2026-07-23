add_test(
  NAME hawaii_pdsim
  COMMAND DynamicSimulation ${CMAKE_CURRENT_SOURCE_DIR}/hawaii.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
