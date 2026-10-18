if(GRIDKIT_ENABLE_ENZYME)
  add_executable(test_opf_components runComponentTests.cpp)
  target_link_libraries(
    test_opf_components
    GridKit::optimal_power_flow_systemmodel
    GridKit::optimal_power_flow_systemmodel_dependency_tracking
    GridKit::testing)
  add_test(NAME OptimalPowerFlowComponentTest COMMAND test_opf_components)
  install(TARGETS test_opf_components)

  add_library(test_opf_coupled_component OBJECT CoupledComponentEnzyme.cpp)
  target_link_libraries(test_opf_coupled_component PRIVATE GridKit::optimal_power_flow_systemmodel
                                                           ClangEnzymeFlags)
  target_compile_options(test_opf_coupled_component PRIVATE -mllvm -enzyme-auto-sparsity=1
                                                            -fno-math-errno)

  add_executable(test_opf_system runSystemTests.cpp)
  target_link_libraries(
    test_opf_system
    test_opf_coupled_component
    GridKit::optimal_power_flow_systemmodel
    GridKit::optimal_power_flow_systemmodel_dependency_tracking
    GridKit::testing)
  add_test(NAME OptimalPowerFlowSystemTest COMMAND test_opf_system)
  install(TARGETS test_opf_system)

  if(GRIDKIT_ENABLE_IPOPT)
    add_executable(test_opf_ipopt runIpoptTests.cpp)
    target_link_libraries(
      test_opf_ipopt
      GridKit::optimal_power_flow_systemmodel
      GridKit::solvers_opt_steady
      GridKit::testing)
    add_test(NAME OptimalPowerFlowIpoptTest COMMAND test_opf_ipopt)
    install(TARGETS test_opf_ipopt)
  endif()
endif()
