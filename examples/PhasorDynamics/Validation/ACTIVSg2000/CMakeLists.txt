gridkit_example_add_file(${CMAKE_SOURCE_DIR}/cases/PhasorDynamics/ACTIVSg2000/ACTIVSg2000.case.json)
gridkit_example_add_file(ACTIVSg2000.solver.json)
add_subdirectory(reference)

add_test(
  NAME ACTIVSg2000_validation
  COMMAND DynamicSimulation ACTIVSg2000.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
