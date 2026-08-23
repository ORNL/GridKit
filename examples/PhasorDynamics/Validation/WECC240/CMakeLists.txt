gridkit_example_add_file(${CMAKE_SOURCE_DIR}/cases/PhasorDynamics/WECC240/WECC240.case.json)
gridkit_example_add_file(WECC240.solver.json)
add_subdirectory(reference)

add_test(
  NAME WECC240_validation
  COMMAND DynamicSimulation WECC240.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
