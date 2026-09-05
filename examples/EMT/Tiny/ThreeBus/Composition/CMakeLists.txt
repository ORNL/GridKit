gridkit_example_add_file(${CMAKE_SOURCE_DIR}/cases/EMT/Composition/TwoSubsystem.case.json)
gridkit_example_add_file(TwoSubsystem.solver.json)
gridkit_example_add_file(${CMAKE_SOURCE_DIR}/cases/EMT/Composition/GovernedMachine.case.json)
gridkit_example_add_file(GovernedMachine.solver.json)

add_test(
  NAME EMTTwoSubsystem
  COMMAND EMTDynamicSimulation TwoSubsystem.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

add_test(
  NAME EMTNestedGovernor
  COMMAND EMTDynamicSimulation GovernedMachine.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
