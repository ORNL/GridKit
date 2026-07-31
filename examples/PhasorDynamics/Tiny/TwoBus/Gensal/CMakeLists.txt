gridkit_example_current_install_path(_install_path)

gridkit_example_add_file(TwoBusGensal.case.json)
gridkit_example_add_file(TwoBusGensal.state.json)
gridkit_example_add_file(TwoBusGensal.solver.json)
gridkit_example_add_file(TwoBusGensal.ref.csv)
add_test(
  NAME TwoBusGensal
  COMMAND DynamicSimulation TwoBusGensal.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

gridkit_example_current_install_path(_install_path)
install(FILES TwoBusGensal.case.json TwoBusGensal.solver.json TwoBusGensal.ref.csv
        DESTINATION ${_install_path})
