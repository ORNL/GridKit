

add_executable(test_phasor_bus TestBus.cpp)
target_link_libraries(test_phasor_bus GRIDKIT::phasor_dynamics_bus)
                                
add_executable(test_phasor_branch TestBus.cpp)
target_link_libraries(test_phasor_branch GRIDKIT::phasor_dynamics_branch
                                         GRIDKIT::phasor_dynamics_bus)
                                
add_test(NAME PhasorDynamicsBusTest COMMAND $<TARGET_FILE:test_phasor_bus>)
add_test(NAME PhasorDynamicsBranchTest COMMAND $<TARGET_FILE:test_phasor_branch>)

install(TARGETS test_phasor_bus test_phasor_branch RUNTIME DESTINATION bin)
