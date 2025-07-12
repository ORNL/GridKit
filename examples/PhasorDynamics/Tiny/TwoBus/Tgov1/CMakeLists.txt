add_executable(TwoBusTgov1 TwoBusTgov1.cpp)
target_link_libraries(TwoBusTgov1 
	GRIDKIT::phasor_dynamics_components
	GRIDKIT::solvers_dyn)
install(TARGETS TwoBusTgov1 RUNTIME DESTINATION bin)

add_test(NAME GenrouTest1_tgov1 COMMAND $<TARGET_FILE:TwoBusTgov1>)

add_executable(TwoBusTgov1Signal TwoBusTgov1Signal.cpp)
target_link_libraries(TwoBusTgov1Signal 
	GRIDKIT::phasor_dynamics_components
	GRIDKIT::phasor_dynamics_signal
	GRIDKIT::solvers_dyn)
install(TARGETS TwoBusTgov1Signal RUNTIME DESTINATION bin)

add_test(NAME TwoBusTgov1SignalTest COMMAND $<TARGET_FILE:TwoBusTgov1Signal>)
