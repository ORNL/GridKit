gridkit_example_current_install_path(_install_path)

add_executable(texas texas.cpp)
target_link_libraries(texas
    GridKit::phasor_dynamics_core
	GridKit::phasor_dynamics_components
	GridKit::phasor_dynamics_signal
	GridKit::solvers_dyn)
install(TARGETS texas RUNTIME DESTINATION ${_install_path})
gridkit_example_add_file(texas.json)

add_test(NAME texas
    COMMAND texas ${CMAKE_CURRENT_BINARY_DIR}/texas.json)
add_test(NAME texas_no_arg
    COMMAND texas
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
