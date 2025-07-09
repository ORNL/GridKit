add_executable(TwoBusTgov1 TwoBusTgov1.cpp)
target_link_libraries(TwoBusTgov1 
	GRIDKIT::phasor_dynamics_bus 
	GRIDKIT::phasor_dynamics_bus_fault 
	GRIDKIT::phasor_dynamics_branch
	GRIDKIT::phasor_dynamics_governortgov1
	GRIDKIT::phasor_dynamics_genrou
	GRIDKIT::phasor_dynamics_gen_classical
	GRIDKIT::phasor_dynamics_load
	GRIDKIT::solvers_dyn)
install(TARGETS TwoBusTgov1 RUNTIME DESTINATION bin)

add_test(NAME GenrouTest1_tgov1 COMMAND $<TARGET_FILE:TwoBusTgov1>)
