gridkit_example_current_install_path(_install_path)

add_executable(hawaiiJson hawaiiJson.cpp)
target_link_libraries(hawaiiJson
    GridKit::phasor_dynamics_core
	GridKit::phasor_dynamics_components
	GridKit::phasor_dynamics_signal
	GridKit::solvers_dyn)
install(TARGETS hawaiiJson RUNTIME DESTINATION ${_install_path})
gridkit_example_add_file(hawaii.json)

add_test(NAME hawaiiJson
    COMMAND hawaiiJson ${CMAKE_CURRENT_BINARY_DIR}/hawaii.json)
add_test(NAME hawaiiJson_no_arg
    COMMAND hawaiiJson
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
