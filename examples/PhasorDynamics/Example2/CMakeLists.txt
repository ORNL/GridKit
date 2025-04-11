add_executable(phasordynamics_example2 example2.cpp)
target_link_libraries(phasordynamics_example2 
	GRIDKIT::phasor_dynamics_bus 
	GRIDKIT::phasor_dynamics_bus_fault 
	GRIDKIT::phasor_dynamics_branch
	GRIDKIT::phasor_dynamics_genrou
	GRIDKIT::phasor_dynamics_load
	SUNDIALS::sunlinsolklu
	SUNDIALS::core
	SUNDIALS::ida
	SUNDIALS::idas
	SUNDIALS::sunmatrixdense)
install(TARGETS phasordynamics_example2 RUNTIME DESTINATION bin)

add_test(NAME GenrouTest COMMAND $<TARGET_FILE:phasordynamics_example2>)
