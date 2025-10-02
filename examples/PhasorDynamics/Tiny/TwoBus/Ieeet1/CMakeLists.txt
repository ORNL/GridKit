add_executable(TwoBusIeeet1 TwoBusIeeet1.cpp)
target_link_libraries(TwoBusIeeet1
	GRIDKIT::phasor_dynamics_components
	GRIDKIT::phasor_dynamics_signal
	GRIDKIT::solvers_dyn)
install(TARGETS TwoBusIeeet1 RUNTIME DESTINATION bin)

add_executable(TwoBusIeeet1Json TwoBusIeeet1Json.cpp)
target_link_libraries(TwoBusIeeet1Json
	GRIDKIT::phasor_dynamics_components
	GRIDKIT::phasor_dynamics_signal
	GRIDKIT::solvers_dyn)
target_include_directories(TwoBusIeeet1Json
    PRIVATE ${CMAKE_SOURCE_DIR}/third-party/nlohmann-json/include)
target_include_directories(TwoBusIeeet1Json
    PRIVATE ${CMAKE_SOURCE_DIR}/third-party/magic-enum/include)
install(TARGETS TwoBusIeeet1Json RUNTIME DESTINATION bin)
gridkit_example_add_file(TwoBusIeeet1.json)

add_test(NAME GenrouTest1_Ieeet1 COMMAND TwoBusIeeet1)
add_test(NAME GenrouTest1_Ieeet1_Json
    COMMAND TwoBusIeeet1Json TwoBusIeeet1.json
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
