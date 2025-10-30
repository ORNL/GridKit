# [[
#  Author(s):
#    - Slaven Peles <peless@ornl.gov>
#]]


add_executable(test_phasor_bus runBusTests.cpp)
target_link_libraries(test_phasor_bus GridKit::phasor_dynamics_bus)

add_executable(test_phasor_branch runBranchTests.cpp)
target_link_libraries(test_phasor_branch GridKit::phasor_dynamics_branch
                                         GridKit::phasor_dynamics_branch_dependency_tracking
                                         GridKit::phasor_dynamics_bus
                                         GridKit::phasor_dynamics_bus_dependency_tracking)

add_executable(test_phasor_load runLoadTests.cpp)
target_link_libraries(test_phasor_load GridKit::phasor_dynamics_load
                                       GridKit::phasor_dynamics_load_dependency_tracking
                                       GridKit::phasor_dynamics_bus
                                       GridKit::phasor_dynamics_bus_dependency_tracking)

add_executable(test_phasor_genrou runGenrouTests.cpp)
target_link_libraries(test_phasor_genrou GridKit::phasor_dynamics_genrou
                                         GridKit::phasor_dynamics_genrou_dependency_tracking
                                         GridKit::phasor_dynamics_bus
                                         GridKit::phasor_dynamics_bus_dependency_tracking)

add_executable(test_phasor_governortgov1 runGovernorTgov1Tests.cpp)
target_link_libraries(test_phasor_governortgov1 GridKit::phasor_dynamics_genrou
                                                GridKit::phasor_dynamics_genrou_dependency_tracking
                                                GridKit::phasor_dynamics_governortgov1
                                                GridKit::phasor_dynamics_governortgov1_dependency_tracking
                                                GridKit::phasor_dynamics_signal
                                                GridKit::phasor_dynamics_bus
                                                GridKit::phasor_dynamics_bus_dependency_tracking)

add_executable(test_phasor_gen_classical runGenClassicalTests.cpp)
target_link_libraries(test_phasor_gen_classical GridKit::phasor_dynamics_gen_classical
                                                GridKit::phasor_dynamics_gen_classical_dependency_tracking
                                                GridKit::phasor_dynamics_bus
                                                GridKit::phasor_dynamics_bus_dependency_tracking)


add_executable(test_phasor_system runSystemTests.cpp)
target_link_libraries(test_phasor_system GridKit::phasor_dynamics_components
                                         GridKit::phasor_dynamics_components_dependency_tracking)

add_executable(test_phasor_system_single_component runSystemSingleComponentTests.cpp)
target_link_libraries(test_phasor_system_single_component GridKit::phasor_dynamics_components
                                                          GridKit::phasor_dynamics_components_dependency_tracking)

add_test(NAME PhasorDynamicsBusTest COMMAND $<TARGET_FILE:test_phasor_bus>)
add_test(NAME PhasorDynamicsBranchTest COMMAND $<TARGET_FILE:test_phasor_branch>)
add_test(NAME PhasorDynamicsGenrouTest COMMAND $<TARGET_FILE:test_phasor_genrou>)
add_test(NAME PhasorDynamicsGovernorTgov1Test COMMAND $<TARGET_FILE:test_phasor_governortgov1>)
add_test(NAME PhasorDynamicsGenClassicalTest COMMAND $<TARGET_FILE:test_phasor_gen_classical>)
add_test(NAME PhasorDynamicsLoadTest COMMAND $<TARGET_FILE:test_phasor_load>)
add_test(NAME PhasorDynamicsSystemTest COMMAND $<TARGET_FILE:test_phasor_system>)
add_test(NAME PhasorDynamicsSystemSingleComponentTest COMMAND $<TARGET_FILE:test_phasor_system_single_component>)

install(TARGETS test_phasor_bus
                test_phasor_branch
                test_phasor_load
                test_phasor_genrou
                test_phasor_governortgov1
                test_phasor_gen_classical
                test_phasor_system 
                test_phasor_system_single_component
        RUNTIME DESTINATION bin)
