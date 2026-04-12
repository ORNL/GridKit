gridkit_example_add_file(newengland.json)
gridkit_example_add_file(newengland.solver.json)

add_test(NAME newengland_pdsim
    COMMAND PDSim newengland.solver.json
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
