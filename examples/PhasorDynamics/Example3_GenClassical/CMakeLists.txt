add_executable(phasordynamics_example3_genclassical example3_genclassical.cpp)
target_link_libraries(phasordynamics_example3_genclassical
	GRIDKIT::phasor_dynamics_bus
	GRIDKIT::phasor_dynamics_bus_fault 
	GRIDKIT::phasor_dynamics_branch
	GRIDKIT::phasor_dynamics_gen_classical
	GRIDKIT::phasor_dynamics_load
	GRIDKIT::solvers_dyn)
install(TARGETS phasordynamics_example3_genclassical RUNTIME DESTINATION bin)

# Not used for tesing for now.
# add_test(NAME GenrouTest3_genclassical COMMAND $<TARGET_FILE:phasordynamics_example3_genclassical>)
