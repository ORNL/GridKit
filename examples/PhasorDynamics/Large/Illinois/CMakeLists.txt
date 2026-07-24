gridkit_example_add_file(illinois.case.json)
gridkit_example_add_file(illinois.solver.json)

add_test(
  NAME illinois_pdsim
  COMMAND DynamicSimulation illinois.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

gridkit_example_current_install_path(_install_path)
install(FILES illinois.case.json illinois.solver.json DESTINATION ${_install_path})
