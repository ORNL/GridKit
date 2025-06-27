add_executable(phasordynamics_example_twobus_tgov1 Example_TwoBus_TGov1.cpp)
target_link_libraries(phasordynamics_example_twobus_tgov1 
	GRIDKIT::phasor_dynamics_bus 
	GRIDKIT::phasor_dynamics_bus_fault 
	GRIDKIT::phasor_dynamics_branch
	GRIDKIT::phasor_dynamics_governortgov1
	GRIDKIT::phasor_dynamics_genrou
	GRIDKIT::phasor_dynamics_load
	GRIDKIT::solvers_dyn)
install(TARGETS phasordynamics_example_twobus_tgov1 RUNTIME DESTINATION bin)

add_test(NAME GenrouTest1_tgov1 COMMAND $<TARGET_FILE:phasordynamics_example_twobus_tgov1>)
