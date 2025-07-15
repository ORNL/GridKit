# [[
#  Author(s):
#    - Slaven Peles <peless@ornl.gov>
#]]

add_executable(ThreeBusClassical ThreeBusClassical.cpp)
target_link_libraries(ThreeBusClassical 
	GRIDKIT::phasor_dynamics_components
	GRIDKIT::solvers_dyn)
install(TARGETS ThreeBusClassical RUNTIME DESTINATION bin)

add_test(NAME ThreeBusClassicalTest COMMAND $<TARGET_FILE:ThreeBusClassical>)
