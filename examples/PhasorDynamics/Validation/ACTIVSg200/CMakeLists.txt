gridkit_example_add_file(${CMAKE_SOURCE_DIR}/cases/PhasorDynamics/ACTIVSg200/ACTIVSg200.case.json)
gridkit_example_add_file(ACTIVSg200.solver.json)
add_subdirectory(reference)

add_test(
  NAME ACTIVSg200_validation
  COMMAND DynamicSimulation ACTIVSg200.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
