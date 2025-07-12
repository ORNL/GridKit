add_executable(TwoBusBasic TwoBusBasic.cpp)
target_link_libraries(TwoBusBasic 
	GRIDKIT::phasor_dynamics_components
	GRIDKIT::solvers_dyn)
install(TARGETS TwoBusBasic RUNTIME DESTINATION bin)

add_test(NAME TwoBusBasicTest COMMAND $<TARGET_FILE:TwoBusBasic>)
