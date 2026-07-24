gridkit_example_add_file(hawaii.case.json)
gridkit_example_add_file(hawaii.solver.json)

add_test(
  NAME hawaii_pdsim
  COMMAND DynamicSimulation hawaii.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

gridkit_example_current_install_path(_install_path)
install(FILES hawaii.case.json hawaii.solver.json DESTINATION ${_install_path})
