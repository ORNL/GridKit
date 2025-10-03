add_executable(TwoBusBasic TwoBusBasic.cpp)
target_link_libraries(TwoBusBasic
	GRIDKIT::phasor_dynamics_components
	GRIDKIT::solvers_dyn)
install(TARGETS TwoBusBasic RUNTIME DESTINATION bin)

add_executable(TwoBusBasicJson TwoBusBasicJson.cpp)
target_link_libraries(TwoBusBasicJson 
	GRIDKIT::phasor_dynamics_components
	GRIDKIT::solvers_dyn)
target_include_directories(TwoBusBasicJson
    PRIVATE ${CMAKE_SOURCE_DIR}/third-party/nlohmann-json/include)
target_include_directories(TwoBusBasicJson
    PRIVATE ${CMAKE_SOURCE_DIR}/third-party/magic-enum/include)
install(TARGETS TwoBusBasicJson RUNTIME DESTINATION bin)
gridkit_example_add_file(TwoBusBasic.json)

add_test(NAME TwoBusBasicTest COMMAND TwoBusBasic)
add_test(NAME TwoBusBasicJsonTest
    COMMAND TwoBusBasicJson TwoBusBasic.json
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
