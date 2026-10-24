add_executable(EconomicDispatch EconomicDispatch.cpp)
target_link_libraries(
  EconomicDispatch
  PUBLIC GridKit::optimal_power_flow_systemmodel GridKit::solvers_opt_steady
         GridKit::utilities_cli_args)
install(TARGETS EconomicDispatch RUNTIME DESTINATION bin)

# Optimal costs [$/h] from PYPOWER 5.1.21
set(_reference_cost_IEEE39 41864.1776)
set(_reference_cost_Hawaii 4151.35488)
set(_reference_cost_ACTIVSg2000 1228892.076)

foreach(case IN ITEMS IEEE39 Hawaii ACTIVSg2000)
  gridkit_example_add_file(${CMAKE_SOURCE_DIR}/cases/OptimalPowerFlow/${case}/${case}.case.json)
  add_test(
    NAME EconomicDispatch_${case}
    COMMAND
      EconomicDispatch --case-file ${case}.case.json --reference-cost ${_reference_cost_${case}}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endforeach()
