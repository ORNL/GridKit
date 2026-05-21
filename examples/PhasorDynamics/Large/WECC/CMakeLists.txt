gridkit_example_add_file(wecc.json)
gridkit_example_add_file(wecc.solver.json)

add_test(NAME wecc_pdsim
    COMMAND DynamicSimulation wecc.solver.json
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
