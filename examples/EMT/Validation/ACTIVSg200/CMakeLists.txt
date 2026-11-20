gridkit_example_add_file(${CMAKE_SOURCE_DIR}/cases/EMT/ACTIVSg200/ACTIVSg200.case.json)
gridkit_example_add_file(${CMAKE_SOURCE_DIR}/cases/EMT/ACTIVSg200/ACTIVSg200.state.json)
gridkit_example_add_file(ACTIVSg200.solver.json)

add_test(
  NAME EMT_ACTIVSg200_validation
  COMMAND TransientSimulation ACTIVSg200.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
