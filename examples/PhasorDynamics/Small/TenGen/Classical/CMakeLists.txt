add_executable(TenGenClassical TenGenClassical.cpp)
target_link_libraries(TenGenClassical
	GridKit::phasor_dynamics_components
	GridKit::solvers_dyn)
install(TARGETS TenGenClassical RUNTIME DESTINATION bin)

# Not used for tesing for now.
# add_test(NAME GenrouTest3_genclassical COMMAND $<TARGET_FILE:TenGenClassical>)
