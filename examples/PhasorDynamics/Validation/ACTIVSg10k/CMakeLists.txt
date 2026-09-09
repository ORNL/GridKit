gridkit_example_add_file(${CMAKE_SOURCE_DIR}/cases/PhasorDynamics/ACTIVSg10k/ACTIVSg10k.case.json)
gridkit_example_add_file(ACTIVSg10k.solver.json)
add_subdirectory(reference)

add_test(
  NAME ACTIVSg10k_validation
  COMMAND DynamicSimulation ACTIVSg10k.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
