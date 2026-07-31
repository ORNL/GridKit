gridkit_example_add_file(wecc.case.json)
gridkit_example_add_file(wecc.state.json)
gridkit_example_add_file(wecc.solver.json)

add_test(
  NAME wecc_pdsim
  COMMAND DynamicSimulation wecc.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

gridkit_example_current_install_path(_install_path)
install(FILES wecc.case.json wecc.solver.json DESTINATION ${_install_path})
