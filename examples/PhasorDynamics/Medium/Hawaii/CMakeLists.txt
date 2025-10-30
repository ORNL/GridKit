gridkit_example_current_install_path(_install_path)

add_executable(hawaiiJson hawaiiJson.cpp)
target_link_libraries(hawaiiJson
	GRIDKIT::phasor_dynamics_components
	GRIDKIT::phasor_dynamics_signal
	GRIDKIT::solvers_dyn)
target_include_directories(hawaiiJson
    PRIVATE ${CMAKE_SOURCE_DIR}/third-party/nlohmann-json/include)
target_include_directories(hawaiiJson
    PRIVATE ${CMAKE_SOURCE_DIR}/third-party/magic-enum/include)
install(TARGETS hawaiiJson RUNTIME DESTINATION ${_install_path})
gridkit_example_add_file(hawaii.json)

add_test(NAME hawaiiJson
    COMMAND hawaiiJson ${CMAKE_CURRENT_BINARY_DIR}/hawaii.json)
add_test(NAME hawaiiJson_no_arg
    COMMAND hawaiiJson
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
