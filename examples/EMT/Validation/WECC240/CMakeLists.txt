gridkit_example_add_file(${CMAKE_SOURCE_DIR}/cases/EMT/WECC240/WECC240.case.json)
gridkit_example_add_file(${CMAKE_SOURCE_DIR}/cases/EMT/WECC240/WECC240.state.json)
gridkit_example_add_file(WECC240.solver.json)

add_test(
  NAME EMT_WECC240_validation
  COMMAND TransientSimulation WECC240.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
