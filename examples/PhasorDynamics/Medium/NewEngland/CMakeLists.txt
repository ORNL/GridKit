gridkit_example_current_install_path(_install_path)

add_executable(newengland newengland.cpp)
target_link_libraries(newengland
    GridKit::phasor_dynamics_core
	GridKit::phasor_dynamics_components
	GridKit::phasor_dynamics_signal
	GridKit::solvers_dyn)
install(TARGETS newengland RUNTIME DESTINATION ${_install_path})
gridkit_example_add_file(newengland.json)

add_test(NAME newengland
    COMMAND newengland ${CMAKE_CURRENT_BINARY_DIR}/newengland.json)
add_test(NAME newengland_no_arg
    COMMAND newengland
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
