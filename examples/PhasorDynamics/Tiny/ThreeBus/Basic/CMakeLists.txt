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
    GridKit::testing)
install(TARGETS ThreeBusBasicJson RUNTIME DESTINATION ${_install_path})
gridkit_example_add_file(ThreeBusBasic.case.json)
gridkit_example_add_file(ThreeBusBasic.ref.csv)
gridkit_example_add_file(ThreeBusBasicBranchTrip.ref.csv)

add_test(NAME ThreeBusBasic COMMAND ThreeBusBasic)
add_test(NAME ThreeBusBasicJson
    COMMAND ThreeBusBasicJson ${CMAKE_CURRENT_BINARY_DIR}/ThreeBusBasic.case.json)
add_test(NAME ThreeBusBasicJson_no_arg
    COMMAND ThreeBusBasicJson
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

gridkit_example_add_file(ThreeBusBasic.solver.json)
gridkit_example_add_file(ThreeBusBasicBranchTrip.solver.json)
add_test(NAME ThreeBusBasic_using_app
    COMMAND PDSim ThreeBusBasic.solver.json
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
add_test(NAME ThreeBusBasicBranchTrip_using_app
    COMMAND PDSim ThreeBusBasicBranchTrip.solver.json
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
