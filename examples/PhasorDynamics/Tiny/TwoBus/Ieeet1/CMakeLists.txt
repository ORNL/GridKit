gridkit_example_add_file(TwoBusIeeet1.case.json)
gridkit_example_add_file(TwoBusIeeet1.state.json)
gridkit_example_add_file(TwoBusIeeet1.solver.json)
gridkit_example_add_file(TwoBusIeeet1.ref.csv)

add_test(
  NAME TwoBusIeeet1
  COMMAND DynamicSimulation TwoBusIeeet1.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

gridkit_example_current_install_path(_install_path)
install(FILES TwoBusIeeet1.case.json TwoBusIeeet1.solver.json TwoBusIeeet1.ref.csv
        DESTINATION ${_install_path})
