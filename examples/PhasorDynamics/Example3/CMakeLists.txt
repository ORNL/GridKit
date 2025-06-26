add_executable(phasordynamics_example3 example3.cpp)
target_link_libraries(phasordynamics_example3
	GRIDKIT::phasor_dynamics_bus
	GRIDKIT::phasor_dynamics_bus_fault 
	GRIDKIT::phasor_dynamics_branch
	GRIDKIT::phasor_dynamics_genrou
	GRIDKIT::phasor_dynamics_load
	GRIDKIT::solvers_dyn)
install(TARGETS phasordynamics_example3 RUNTIME DESTINATION bin)

# Not used for tesing for now.
# add_test(NAME GenrouTest3 COMMAND $<TARGET_FILE:phasordynamics_example3>)
