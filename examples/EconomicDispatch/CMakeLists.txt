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
# gridkit_economic_dispatch_example(<case>)
#
# Solve the economic dispatch of `cases/PhasorDynamics/<case>` with the limits
# and costs of its MATPOWER case, and write the solution state.
#
function(gridkit_economic_dispatch_example case)
  gridkit_example_add_file(${CMAKE_SOURCE_DIR}/cases/PhasorDynamics/${case}/${case}.case.json)
  gridkit_example_add_file(${CMAKE_SOURCE_DIR}/cases/PhasorDynamics/${case}/${case}.m)
  gridkit_example_add_file(${case}.solver.json)

  add_test(
    NAME EconomicDispatch_${case}
    COMMAND EconomicDispatch ${case}.solver.json
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endfunction()

add_subdirectory(IEEE39)
add_subdirectory(Hawaii)
add_subdirectory(ACTIVSg2000)
