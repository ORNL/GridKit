gridkit_example_add_file(IEEE39.A.solver.json)
gridkit_example_add_file(IEEE39.B.solver.json)
gridkit_example_add_file(${CMAKE_SOURCE_DIR}/cases/PhasorDynamics/WECC240/WECC240.case.json)
gridkit_example_add_file(${CMAKE_SOURCE_DIR}/cases/PhasorDynamics/WECC240/WECC240.state.json)
gridkit_example_add_file(WECC240.monolithic.solver.json)
set(wecc_regions
    SoCal
    NorCal
    Southwest
    Northwest)
foreach(region IN LISTS wecc_regions)
  gridkit_example_add_file(WECC240.${region}.solver.json)
endforeach()
add_subdirectory(partition)
add_subdirectory(reference)

set(simulate $<TARGET_FILE:DynamicSimulation>)
add_test(
  NAME IEEE39_cosimulation
  COMMAND sh -c "${simulate} IEEE39.A.solver.json & ${simulate} IEEE39.B.solver.json && wait $!"
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

add_test(
  NAME WECC240_cosimulation
  COMMAND
    sh -c [=[
    pids=""
    for region do
      "$0" "WECC240.$region.solver.json" &
      pids="$pids $!"
    done
    status=0
    for pid in $pids; do wait "$pid" || status=1; done
    exit "$status"
  ]=] ${simulate} ${wecc_regions}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
