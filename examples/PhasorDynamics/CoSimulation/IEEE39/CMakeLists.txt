foreach(_side Grid Generator)
  gridkit_example_add_file(
    ${CMAKE_SOURCE_DIR}/cases/PhasorDynamics/IEEE39/IEEE39${_side}.case.json)
  gridkit_example_add_file(IEEE39${_side}.solver.json)
  gridkit_example_add_file(IEEE39${_side}.contingency.solver.json)
endforeach()
