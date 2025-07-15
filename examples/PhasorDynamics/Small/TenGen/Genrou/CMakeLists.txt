add_executable(TenGenGenrou TenGenGenrou.cpp)
target_link_libraries(TenGenGenrou
	GRIDKIT::phasor_dynamics_components
	GRIDKIT::solvers_dyn)
install(TARGETS TenGenGenrou RUNTIME DESTINATION bin)

# Not used for tesing for now.
# add_test(NAME GenrouTest3 COMMAND $<TARGET_FILE:TenGenGenrou>)
