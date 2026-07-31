gridkit_example_add_file(newengland.case.json)
gridkit_example_add_file(newengland.state.json)
gridkit_example_add_file(newengland.solver.json)

add_test(
  NAME newengland_pdsim
  COMMAND DynamicSimulation newengland.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

add_test(
  NAME newengland_ca
  COMMAND ContingencyAnalysis newengland.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

gridkit_example_current_install_path(_install_path)
install(FILES newengland.case.json newengland.solver.json DESTINATION ${_install_path})
