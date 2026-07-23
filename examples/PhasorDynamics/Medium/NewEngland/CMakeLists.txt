add_test(
  NAME newengland_pdsim
  COMMAND DynamicSimulation ${CMAKE_CURRENT_SOURCE_DIR}/newengland.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

add_test(
  NAME newengland_ca
  COMMAND ContingencyAnalysis ${CMAKE_CURRENT_SOURCE_DIR}/newengland.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
