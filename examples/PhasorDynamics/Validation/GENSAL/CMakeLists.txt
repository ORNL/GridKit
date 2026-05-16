gridkit_example_add_file(GENSAL.case.json)
gridkit_example_add_file(GENSAL.solver.json)
gridkit_example_add_file(GENSAL.ref.csv)

add_test(NAME GENSAL_validation
    COMMAND PDSim GENSAL.solver.json
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(GENSAL_validation PROPERTIES
    LABELS "validation"
    TIMEOUT 60)
