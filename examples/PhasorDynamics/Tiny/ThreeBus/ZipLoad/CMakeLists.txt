gridkit_example_add_file(ThreeBusZipLoad.case.json)
gridkit_example_add_file(ThreeBusZipLoad.solver.json)
gridkit_example_add_file(ThreeBusZipLoad.ref.csv)

add_test(
  NAME ThreeBusZipLoad
  COMMAND DynamicSimulation ThreeBusZipLoad.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

gridkit_example_current_install_path(_install_path)
install(FILES ThreeBusZipLoad.case.json ThreeBusZipLoad.solver.json ThreeBusZipLoad.ref.csv
        DESTINATION ${_install_path})
