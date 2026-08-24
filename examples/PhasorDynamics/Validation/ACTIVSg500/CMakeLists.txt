gridkit_example_add_file(${CMAKE_SOURCE_DIR}/cases/PhasorDynamics/ACTIVSg500/ACTIVSg500.case.json)
gridkit_example_add_file(ACTIVSg500.solver.json)
add_subdirectory(reference)

add_test(
  NAME ACTIVSg500_validation
  COMMAND DynamicSimulation ACTIVSg500.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
