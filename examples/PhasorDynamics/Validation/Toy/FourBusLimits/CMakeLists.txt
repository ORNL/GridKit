gridkit_example_add_file(${CMAKE_SOURCE_DIR}/cases/PhasorDynamics/Toy/FourBusLimits.case.json)
gridkit_example_add_file(${CMAKE_SOURCE_DIR}/cases/PhasorDynamics/Toy/FourBusLimits.state.json)
gridkit_example_add_file(FourBusLimits.solver.json)
add_subdirectory(reference)

add_test(
  NAME FourBusLimits_validation
  COMMAND DynamicSimulation FourBusLimits.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
