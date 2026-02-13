gridkit_example_current_install_path(_install_path)

add_executable(ThreeBusBasic ThreeBusBasic.cpp)
target_link_libraries(ThreeBusBasic
	GridKit::phasor_dynamics_components
	GridKit::solvers_dyn)
install(TARGETS ThreeBusBasic RUNTIME DESTINATION ${_install_path})

add_executable(ThreeBusBasicJson ThreeBusBasicJson.cpp)
target_link_libraries(ThreeBusBasicJson
	GridKit::phasor_dynamics_components
	GridKit::solvers_dyn
    GridKit::Utilities)
install(TARGETS ThreeBusBasicJson RUNTIME DESTINATION ${_install_path})
gridkit_example_add_file(ThreeBusBasic.json)
gridkit_example_add_file(ThreeBusBasic.ref.csv)

add_test(NAME ThreeBusBasic COMMAND ThreeBusBasic)
add_test(NAME ThreeBusBasicJson
    COMMAND ThreeBusBasicJson ${CMAKE_CURRENT_BINARY_DIR}/ThreeBusBasic.json)
add_test(NAME ThreeBusBasicJson_no_arg
    COMMAND ThreeBusBasicJson
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
