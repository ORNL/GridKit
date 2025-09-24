add_executable(TwoBusTgov1 TwoBusTgov1.cpp)
target_link_libraries(TwoBusTgov1 
	GRIDKIT::phasor_dynamics_components
	GRIDKIT::phasor_dynamics_signal
	GRIDKIT::solvers_dyn)
install(TARGETS TwoBusTgov1 RUNTIME DESTINATION bin)

add_executable(TwoBusTgov1Json TwoBusTgov1Json.cpp)
target_link_libraries(TwoBusTgov1Json 
	GRIDKIT::phasor_dynamics_components
	GRIDKIT::phasor_dynamics_signal
	GRIDKIT::solvers_dyn)
target_include_directories(TwoBusTgov1Json PRIVATE ${CMAKE_SOURCE_DIR}/third-party/nlohmann-json/include)
target_include_directories(TwoBusTgov1Json PRIVATE ${CMAKE_SOURCE_DIR}/third-party/magic-enum/include)
install(TARGETS TwoBusTgov1Json RUNTIME DESTINATION bin)

add_test(NAME GenrouTest1_tgov1 COMMAND $<TARGET_FILE:TwoBusTgov1>)
add_test(NAME GenrouTest1_tgov1_json COMMAND $<TARGET_FILE:TwoBusTgov1Json>)
