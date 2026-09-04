gridkit_example_add_file(${CMAKE_SOURCE_DIR}/cases/EMT/OverheadLine/OverheadLine.case.json)
gridkit_example_add_file(OverheadLine.solver.json)

add_test(
  NAME EMTOverheadLine
  COMMAND EMTDynamicSimulation OverheadLine.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
