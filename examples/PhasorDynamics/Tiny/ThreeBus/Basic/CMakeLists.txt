add_executable(ThreeBusBasic ThreeBusBasic.cpp)
target_link_libraries(ThreeBusBasic 
	GRIDKIT::phasor_dynamics_bus 
	GRIDKIT::phasor_dynamics_bus_fault 
	GRIDKIT::phasor_dynamics_branch
	GRIDKIT::phasor_dynamics_genrou
	GRIDKIT::phasor_dynamics_load
	GRIDKIT::solvers_dyn)
install(TARGETS ThreeBusBasic RUNTIME DESTINATION bin)

add_test(NAME GenrouTest2 COMMAND $<TARGET_FILE:ThreeBusBasic>)
