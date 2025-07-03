add_executable(TenGenGenrou TenGenGenrou.cpp)
target_link_libraries(TenGenGenrou
	GRIDKIT::phasor_dynamics_bus 
	GRIDKIT::phasor_dynamics_bus_fault 
	GRIDKIT::phasor_dynamics_branch
	GRIDKIT::phasor_dynamics_genrou
	GRIDKIT::phasor_dynamics_governortgov1
	GRIDKIT::phasor_dynamics_load
	GRIDKIT::solvers_dyn)
install(TARGETS TenGenGenrou RUNTIME DESTINATION bin)

# Not used for tesing for now.
# add_test(NAME GenrouTest3 COMMAND $<TARGET_FILE:TenGenGenrou>)
