gridkit_example_current_install_path(_install_path)

gridkit_example_add_file(Gensal.case.json)
gridkit_example_add_file(Gensal.solver.json)
gridkit_example_add_file(Gensal.ref.csv)
add_test(NAME TwoBusGensal
    COMMAND PDSim Gensal.solver.json
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
