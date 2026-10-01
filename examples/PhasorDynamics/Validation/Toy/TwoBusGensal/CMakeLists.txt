gridkit_example_add_file(${CMAKE_SOURCE_DIR}/cases/PhasorDynamics/Toy/TwoBusGensal.case.json)
gridkit_example_add_file(TwoBusGensal.solver.json)
add_subdirectory(reference)

add_test(
  NAME TwoBusGensal_validation
  COMMAND DynamicSimulation TwoBusGensal.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
