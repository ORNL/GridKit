gridkit_example_add_file(${CMAKE_SOURCE_DIR}/cases/PhasorDynamics/Toy/ThreeBusZipLoad.case.json)
gridkit_example_add_file(ThreeBusZipLoad.solver.json)
add_subdirectory(reference)

add_test(
  NAME ThreeBusZipLoad_validation
  COMMAND DynamicSimulation ThreeBusZipLoad.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
