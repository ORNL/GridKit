add_executable(TwoBusTgov1 TwoBusTgov1.cpp)
target_link_libraries(TwoBusTgov1 
	GRIDKIT::phasor_dynamics_components
	GRIDKIT::phasor_dynamics_signal
	GRIDKIT::solvers_dyn)
install(TARGETS TwoBusTgov1 RUNTIME DESTINATION bin)

add_test(NAME GenrouTest1_tgov1 COMMAND $<TARGET_FILE:TwoBusTgov1>)

