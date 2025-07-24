add_executable(TwoBusIeeet1 TwoBusIeeet1.cpp)
target_link_libraries(TwoBusIeeet1 
	GRIDKIT::phasor_dynamics_components
	GRIDKIT::phasor_dynamics_signal
	GRIDKIT::solvers_dyn)
install(TARGETS TwoBusIeeet1 RUNTIME DESTINATION bin)

add_test(NAME GenrouTest1_Ieeet1 COMMAND $<TARGET_FILE:TwoBusIeeet1>)
