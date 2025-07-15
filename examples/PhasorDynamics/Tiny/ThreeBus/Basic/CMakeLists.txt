add_executable(ThreeBusBasic ThreeBusBasic.cpp)
target_link_libraries(ThreeBusBasic
	GRIDKIT::phasor_dynamic_components
	GRIDKIT::solvers_dyn)
install(TARGETS ThreeBusBasic RUNTIME DESTINATION bin)

add_test(NAME GenrouTest2 COMMAND $<TARGET_FILE:ThreeBusBasic>)
