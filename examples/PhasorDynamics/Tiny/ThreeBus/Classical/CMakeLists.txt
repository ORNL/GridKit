# [[
#  Author(s):
#    - Slaven Peles <peless@ornl.gov>
#]]

add_executable(ThreeBusClassical ThreeBusClassical.cpp)
target_link_libraries(ThreeBusClassical 
	GRIDKIT::phasor_dynamics_bus 
	GRIDKIT::phasor_dynamics_bus_fault 
	GRIDKIT::phasor_dynamics_branch
	GRIDKIT::phasor_dynamics_genrou
	GRIDKIT::phasor_dynamics_gen_classical
	GRIDKIT::phasor_dynamics_load
	GRIDKIT::solvers_dyn)
install(TARGETS ThreeBusClassical RUNTIME DESTINATION bin)

add_test(NAME ThreeBusClassicalTest COMMAND $<TARGET_FILE:ThreeBusClassical>)
