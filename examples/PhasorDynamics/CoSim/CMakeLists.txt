gridkit_example_add_file(IEEE39.A.solver.json)
gridkit_example_add_file(IEEE39.B.solver.json)
add_subdirectory(partition)
add_subdirectory(reference)

set(simulate $<TARGET_FILE:DynamicSimulation>)
add_test(
  NAME IEEE39_cosimulation
  COMMAND sh -c "${simulate} IEEE39.A.solver.json & ${simulate} IEEE39.B.solver.json && wait $!"
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
