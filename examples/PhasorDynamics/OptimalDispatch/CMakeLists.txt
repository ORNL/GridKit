#
# gridkit_optimal_dispatch_example(<case> <residual-tolerance> <drift-tolerance>)
#
# Solve the optimal dispatch of `cases/PhasorDynamics/<case>` with the limits
# and costs of its MATPOWER case, and check that PhasorDynamics starts in
# steady state from the solution. The tolerances bound the residual at t = 0
# and the relative drift over 1 s.
#
function(
  gridkit_optimal_dispatch_example
  case
  residual_tol
  drift_tol)
  gridkit_example_add_file(${CMAKE_SOURCE_DIR}/cases/PhasorDynamics/${case}/${case}.case.json)
  gridkit_example_add_file(${CMAKE_SOURCE_DIR}/cases/PhasorDynamics/${case}/${case}.m)
  gridkit_example_add_file(${case}.solver.json)

  add_test(
    NAME OptimalDispatch_${case}
    COMMAND OptimalDispatch ${case}.solver.json
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties(OptimalDispatch_${case} PROPERTIES FIXTURES_SETUP ${case}_dispatch_state)

  add_test(
    NAME OptimalDispatch_${case}_steady_start
    COMMAND test_opf_steady_start ${case}.case.json ${case}.state.json ${residual_tol} ${drift_tol}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties(OptimalDispatch_${case}_steady_start PROPERTIES FIXTURES_REQUIRED
                                                                       ${case}_dispatch_state)
endfunction()

add_subdirectory(IEEE39)
add_subdirectory(Hawaii)
add_subdirectory(ACTIVSg2000)
