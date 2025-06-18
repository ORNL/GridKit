add_executable(phasordynamics_example3 example3.cpp)
target_link_libraries(phasordynamics_example3
	GRIDKIT::phasor_dynamics_bus 
	GRIDKIT::phasor_dynamics_bus_fault 
	GRIDKIT::phasor_dynamics_branch
	GRIDKIT::phasor_dynamics_genrou
	GRIDKIT::phasor_dynamics_turbinegov
	SUNDIALS::sunlinsolklu
	SUNDIALS::core
	SUNDIALS::ida
	SUNDIALS::idas
	SUNDIALS::sunmatrixdense)
install(TARGETS phasordynamics_example3 RUNTIME DESTINATION bin)

add_test(NAME GenrouTest COMMAND $<TARGET_FILE:phasordynamics_example3>)
