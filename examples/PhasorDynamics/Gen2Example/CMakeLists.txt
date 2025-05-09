add_executable(gen2_example example.cpp)
target_link_libraries(gen2_example 
	GRIDKIT::phasor_dynamics_bus 
	GRIDKIT::phasor_dynamics_bus_fault 
	GRIDKIT::phasor_dynamics_branch
	GRIDKIT::phasor_dynamics_classical_gen
	GRIDKIT::solvers_dyn)
install(TARGETS gen2_example RUNTIME DESTINATION bin)

add_test(NAME GenClassicalTest1 COMMAND $<TARGET_FILE:gen2_example>)

