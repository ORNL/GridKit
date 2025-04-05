add_executable(phasordynamics_example1 example1.cpp)
target_link_libraries(phasordynamics_example1 
	GRIDKIT::phasor_dynamics_bus 
	GRIDKIT::phasor_dynamics_branch
	GRIDKIT::phasor_dynamics_genrou
	SUNDIALS::sunlinsolklu
	SUNDIALS::core
	SUNDIALS::ida
	SUNDIALS::idas
	SUNDIALS::sunmatrixdense)
install(TARGETS phasordynamics_example1 RUNTIME DESTINATION bin)

add_test(NAME GenrouTest COMMAND $<TARGET_FILE:phasordynamics_example1>)
