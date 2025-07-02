add_executable(TwoBusBasic TwoBusBasic.cpp)
target_link_libraries(TwoBusBasic 
	GRIDKIT::phasor_dynamics_bus 
	GRIDKIT::phasor_dynamics_bus_fault 
	GRIDKIT::phasor_dynamics_branch
	GRIDKIT::phasor_dynamics_genrou
	GRIDKIT::phasor_dynamics_gen_classical
	GRIDKIT::phasor_dynamics_load
	GRIDKIT::solvers_dyn)
install(TARGETS TwoBusBasic RUNTIME DESTINATION bin)

add_test(NAME GenrouTest1 COMMAND $<TARGET_FILE:TwoBusBasic>)
