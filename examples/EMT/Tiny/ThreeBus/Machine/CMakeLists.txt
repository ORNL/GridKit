gridkit_example_add_file(Machine.case.json)
gridkit_example_add_file(Machine.state.json)
gridkit_example_add_file(Machine.solver.json)
gridkit_example_add_file(Machine.ref.csv)

add_test(
  NAME EMTMachine
  COMMAND EMTDynamicSimulation Machine.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

gridkit_example_current_install_path(_install_path)
install(
  FILES Machine.case.json
        Machine.state.json
        Machine.solver.json
        Machine.ref.csv
  DESTINATION ${_install_path})
