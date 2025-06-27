add_executable(phasordynamics_example_tengen_genrou Example_TenGen_Genrou.cpp)
target_link_libraries(phasordynamics_example_tengen_genrou
	GRIDKIT::phasor_dynamics_bus 
	GRIDKIT::phasor_dynamics_bus_fault 
	GRIDKIT::phasor_dynamics_branch
	GRIDKIT::phasor_dynamics_genrou
	GRIDKIT::phasor_dynamics_load
	GRIDKIT::solvers_dyn)
install(TARGETS phasordynamics_example_tengen_genrou RUNTIME DESTINATION bin)

# Not used for tesing for now.
# add_test(NAME GenrouTest3 COMMAND $<TARGET_FILE:phasordynamics_example_tengen_genrou>)
