add_executable(test_opf_components runComponentTests.cpp)
target_link_libraries(
  test_opf_components
  GridKit::optimal_power_flow_components
  GridKit::optimal_power_flow_components_dependency_tracking
  GridKit::phasor_dynamics_branch
  GridKit::phasor_dynamics_bus
  GridKit::testing)
add_test(NAME OptimalPowerFlowComponentTest COMMAND test_opf_components)
install(TARGETS test_opf_components)

if(GRIDKIT_ENABLE_IPOPT)
  add_executable(test_opf_system runSystemTests.cpp)
  target_link_libraries(
    test_opf_system
    GridKit::optimal_power_flow_systemmodel
    GridKit::optimal_power_flow_systemmodel_dependency_tracking
    GridKit::solvers_opt_problem
    GridKit::testing)
  add_test(NAME OptimalPowerFlowSystemTest COMMAND test_opf_system)
  install(TARGETS test_opf_system)
endif()
