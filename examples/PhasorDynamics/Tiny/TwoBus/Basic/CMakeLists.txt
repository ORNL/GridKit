gridkit_example_current_install_path(_install_path)

add_executable(TwoBusBasic TwoBusBasic.cpp)
target_link_libraries(TwoBusBasic
	GridKit::phasor_dynamics_components
	GridKit::solvers_dyn)
install(TARGETS TwoBusBasic RUNTIME DESTINATION ${_install_path})

add_executable(TwoBusBasicJson TwoBusBasicJson.cpp)
target_link_libraries(TwoBusBasicJson
	GridKit::phasor_dynamics_components
	GridKit::solvers_dyn)
install(TARGETS TwoBusBasicJson RUNTIME DESTINATION ${_install_path})
gridkit_example_add_file(TwoBusBasic.case.json)

add_test(NAME TwoBusBasic COMMAND TwoBusBasic)
add_test(NAME TwoBusBasicJson
    COMMAND TwoBusBasicJson ${CMAKE_CURRENT_BINARY_DIR}/TwoBusBasic.case.json)
add_test(NAME TwoBusBasicJson_no_arg
    COMMAND TwoBusBasicJson
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

gridkit_example_add_file(TwoBusBasic.solver.json)
gridkit_example_add_file(TwoBusBasic.ref.csv)
add_test(NAME TwoBusBasic_using_app
    COMMAND PDSim TwoBusBasic.solver.json
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
