gridkit_example_add_file(ThreeBusBasic.case.json)
gridkit_example_add_file(ThreeBusBasic.ref.csv)
gridkit_example_add_file(ThreeBusBasic.state.json)
gridkit_example_add_file(ThreeBusBasic.solver.json)

add_test(
  NAME ThreeBusBasic
  COMMAND DynamicSimulation ThreeBusBasic.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

gridkit_example_current_install_path(_install_path)
install(FILES ThreeBusBasic.case.json ThreeBusBasic.solver.json ThreeBusBasic.ref.csv
        DESTINATION ${_install_path})
