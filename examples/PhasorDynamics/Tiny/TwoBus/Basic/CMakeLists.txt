gridkit_example_current_install_path(_install_path)

add_executable(TwoBusBasic TwoBusBasic.cpp)
target_link_libraries(TwoBusBasic
	GRIDKIT::phasor_dynamics_components
	GRIDKIT::solvers_dyn)
install(TARGETS TwoBusBasic RUNTIME DESTINATION ${_install_path})

add_executable(TwoBusBasicJson TwoBusBasicJson.cpp)
target_link_libraries(TwoBusBasicJson
	GRIDKIT::phasor_dynamics_components
	GRIDKIT::solvers_dyn)
install(TARGETS TwoBusBasicJson RUNTIME DESTINATION ${_install_path})
gridkit_example_add_file(TwoBusBasic.json)

add_test(NAME TwoBusBasic COMMAND TwoBusBasic)
add_test(NAME TwoBusBasicJson
    COMMAND TwoBusBasicJson ${CMAKE_CURRENT_BINARY_DIR}/TwoBusBasic.json)
add_test(NAME TwoBusBasicJson_no_arg
    COMMAND TwoBusBasicJson
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
