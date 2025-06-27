add_executable(phasordynamics_example_threebus_basic Example_ThreeBus_Basic.cpp)
target_link_libraries(phasordynamics_example_threebus_basic 
	GRIDKIT::phasor_dynamics_bus 
	GRIDKIT::phasor_dynamics_bus_fault 
	GRIDKIT::phasor_dynamics_branch
	GRIDKIT::phasor_dynamics_genrou
	GRIDKIT::phasor_dynamics_load
	GRIDKIT::solvers_dyn)
install(TARGETS phasordynamics_example_threebus_basic RUNTIME DESTINATION bin)

add_test(NAME GenrouTest2 COMMAND $<TARGET_FILE:phasordynamics_example_threebus_basic>)
