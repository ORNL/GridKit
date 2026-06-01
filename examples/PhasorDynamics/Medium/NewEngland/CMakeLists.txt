gridkit_example_add_file(newengland.json)
gridkit_example_add_file(newengland.solver.json)

add_test(
  NAME newengland_pdsim
  COMMAND DynamicSimulation newengland.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

add_test(
  NAME newengland_ca
  COMMAND ContingencyAnalysis newengland.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(newengland_ca PROPERTIES WILL_FAIL TRUE)
