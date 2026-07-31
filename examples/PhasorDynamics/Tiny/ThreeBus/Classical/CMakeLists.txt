# [[
#  Author(s):
#    - Slaven Peles <peless@ornl.gov>
#]]

gridkit_example_add_file(ThreeBusClassical.case.json)
gridkit_example_add_file(ThreeBusClassical.state.json)
gridkit_example_add_file(ThreeBusClassical.solver.json)
gridkit_example_add_file(ThreeBusClassical.ref.csv)

add_test(
  NAME ThreeBusClassical
  COMMAND DynamicSimulation ThreeBusClassical.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

gridkit_example_current_install_path(_install_path)
install(FILES ThreeBusClassical.case.json ThreeBusClassical.solver.json ThreeBusClassical.ref.csv
        DESTINATION ${_install_path})
