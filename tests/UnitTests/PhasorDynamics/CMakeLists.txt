# [[
#  Author(s):
#    - Slaven Peles <peless@ornl.gov>
#]]


add_executable(test_phasor_bus runBusTests.cpp)
target_link_libraries(test_phasor_bus GRIDKIT::phasor_dynamics_bus)
                                
add_executable(test_phasor_branch runBranchTests.cpp)
target_link_libraries(test_phasor_branch GRIDKIT::phasor_dynamics_branch
                                         GRIDKIT::phasor_dynamics_bus)
                                
add_executable(test_phasor_load runLoadTests.cpp)
target_link_libraries(test_phasor_load GRIDKIT::phasor_dynamics_load
                                       GRIDKIT::phasor_dynamics_bus)
                                
add_executable(test_phasor_system runSystemTests.cpp)
target_link_libraries(test_phasor_system GRIDKIT::phasor_dynamics_load
                                         GRIDKIT::phasor_dynamics_branch
                                         GRIDKIT::phasor_dynamics_bus)

add_test(NAME PhasorDynamicsBusTest COMMAND $<TARGET_FILE:test_phasor_bus>)
add_test(NAME PhasorDynamicsBranchTest COMMAND $<TARGET_FILE:test_phasor_branch>)
add_test(NAME PhasorDynamicsLoadTest COMMAND $<TARGET_FILE:test_phasor_load>)
add_test(NAME PhasorDynamicsSystemTest COMMAND $<TARGET_FILE:test_phasor_system>)

install(TARGETS test_phasor_bus test_phasor_branch test_phasor_load test_phasor_system RUNTIME DESTINATION bin)
