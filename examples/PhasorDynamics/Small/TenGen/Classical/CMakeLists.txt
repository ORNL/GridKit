add_executable(TenGenClassical TenGenClassical.cpp)
target_link_libraries(TenGenClassical
	GRIDKIT::phasor_dynamics_bus 
	GRIDKIT::phasor_dynamics_bus_fault 
	GRIDKIT::phasor_dynamics_branch
	GRIDKIT::phasor_dynamics_gen_classical
	GRIDKIT::phasor_dynamics_governortgov1
	GRIDKIT::phasor_dynamics_load
	GRIDKIT::solvers_dyn)
install(TARGETS TenGenClassical RUNTIME DESTINATION bin)

# Not used for tesing for now.
# add_test(NAME GenrouTest3_genclassical COMMAND $<TARGET_FILE:TenGenClassical>)
