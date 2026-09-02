gridkit_example_add_file(OverheadLine.case.json)
gridkit_example_add_file(OverheadLine.solver.json)
gridkit_example_add_file(OverheadLine.ref.csv)

add_test(
  NAME EMTOverheadLine
  COMMAND EMTDynamicSimulation OverheadLine.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

gridkit_example_current_install_path(_install_path)
install(FILES OverheadLine.case.json OverheadLine.solver.json OverheadLine.ref.csv
        DESTINATION ${_install_path})
