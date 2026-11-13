gridkit_example_add_file(${CMAKE_SOURCE_DIR}/cases/EMT/Toy/ThreeBusBasic.case.json)
gridkit_example_add_file(ThreeBusBasic.solver.json)

add_test(
  NAME EMT_ThreeBusBasic_line_energization
  COMMAND TransientSimulation ThreeBusBasic.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
