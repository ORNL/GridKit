gridkit_example_current_install_path(_install_path)

gridkit_example_add_file(overhead.line.json)
gridkit_example_add_file(overhead.solver.json)
gridkit_example_add_file(yc.solver.json)
gridkit_example_add_file(hmin.solver.json)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/output)
install(DIRECTORY DESTINATION ${_install_path}/output)

add_test(NAME FrequencyResponseOverhead COMMAND FrequencyResponse overhead.solver.json)
