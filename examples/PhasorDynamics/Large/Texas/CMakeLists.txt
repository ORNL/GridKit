gridkit_example_add_file(texas.json)
gridkit_example_add_file(texas.solver.json)

add_test(NAME texas_pdsim
    COMMAND DynamicSimulation texas.solver.json
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
