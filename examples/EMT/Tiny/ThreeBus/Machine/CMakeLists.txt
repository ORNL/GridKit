gridkit_example_add_file(${CMAKE_SOURCE_DIR}/cases/EMT/Machine/Machine.case.json)
gridkit_example_add_file(Machine.state.json)
gridkit_example_add_file(Machine.solver.json)

add_test(
  NAME EMTMachine
  COMMAND EMTDynamicSimulation Machine.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
