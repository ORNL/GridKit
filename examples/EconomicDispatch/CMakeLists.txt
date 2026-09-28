add_executable(EconomicDispatch EconomicDispatch.cpp)
target_link_libraries(
  EconomicDispatch
  PUBLIC GridKit::optimal_power_flow_systemmodel
         GridKit::phasor_dynamics_systemmodel
         GridKit::model_state_data
         GridKit::solvers_opt_steady)
target_include_directories(
  EconomicDispatch
  PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/nlohmann-json/include)
install(
  TARGETS EconomicDispatch
  EXPORT gridkit-targets
  RUNTIME)

#
# gridkit_economic_dispatch_example(<case> <residual-tolerance> <drift-tolerance>)
#
# Solve the economic dispatch of `cases/PhasorDynamics/<case>` with the limits
# and costs of its MATPOWER case, and check that PhasorDynamics starts in
# steady state from the solution. The tolerances bound the residual at t = 0
# and the relative drift over 1 s.
#
function(
  gridkit_economic_dispatch_example
  case
  residual_tol
  drift_tol)
  gridkit_example_add_file(${CMAKE_SOURCE_DIR}/cases/PhasorDynamics/${case}/${case}.case.json)
  gridkit_example_add_file(${CMAKE_SOURCE_DIR}/cases/PhasorDynamics/${case}/${case}.m)
  gridkit_example_add_file(${case}.solver.json)

  add_test(
    NAME EconomicDispatch_${case}
    COMMAND EconomicDispatch ${case}.solver.json
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties(EconomicDispatch_${case} PROPERTIES FIXTURES_SETUP ${case}_dispatch_state)

  add_test(
    NAME EconomicDispatch_${case}_steady_start
    COMMAND test_opf_steady_start ${case}.case.json ${case}.state.json ${residual_tol} ${drift_tol}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties(EconomicDispatch_${case}_steady_start PROPERTIES FIXTURES_REQUIRED
                                                                        ${case}_dispatch_state)
endfunction()

add_subdirectory(IEEE39)
add_subdirectory(Hawaii)
add_subdirectory(ACTIVSg2000)
