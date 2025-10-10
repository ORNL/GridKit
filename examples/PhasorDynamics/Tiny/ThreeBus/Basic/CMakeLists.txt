gridkit_example_current_install_path(_install_path)

add_executable(ThreeBusBasic ThreeBusBasic.cpp)
target_link_libraries(ThreeBusBasic
	GRIDKIT::phasor_dynamics_components
	GRIDKIT::solvers_dyn)
install(TARGETS ThreeBusBasic RUNTIME DESTINATION ${_install_path})

add_executable(ThreeBusBasicJson ThreeBusBasicJson.cpp)
target_link_libraries(ThreeBusBasicJson
	GRIDKIT::phasor_dynamics_components
	GRIDKIT::solvers_dyn)
target_include_directories(ThreeBusBasicJson
    PRIVATE ${CMAKE_SOURCE_DIR}/third-party/nlohmann-json/include)
target_include_directories(ThreeBusBasicJson
    PRIVATE ${CMAKE_SOURCE_DIR}/third-party/magic-enum/include)
install(TARGETS ThreeBusBasicJson RUNTIME DESTINATION ${_install_path})
gridkit_example_add_file(ThreeBusBasic.json)

add_test(NAME GenrouTest2 COMMAND ThreeBusBasic)
add_test(NAME GenrouTest2Json
    COMMAND ThreeBusBasicJson ${CMAKE_CURRENT_BINARY_DIR}/ThreeBusBasic.json)
add_test(NAME GenrouTest2Json_no_arg
    COMMAND ThreeBusBasicJson
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
