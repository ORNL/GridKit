gridkit_example_add_file(texas.case.json)
gridkit_example_add_file(texas.solver.json)

add_test(
  NAME texas_pdsim
  COMMAND DynamicSimulation texas.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

gridkit_example_current_install_path(_install_path)
install(FILES texas.case.json texas.solver.json DESTINATION ${_install_path})
