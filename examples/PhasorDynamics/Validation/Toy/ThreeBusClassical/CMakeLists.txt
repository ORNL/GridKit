# [[
#  Author(s):
#    - Slaven Peles <peless@ornl.gov>
#]]

gridkit_example_add_file(${CMAKE_SOURCE_DIR}/cases/PhasorDynamics/Toy/ThreeBusClassical.case.json)
gridkit_example_add_file(ThreeBusClassical.solver.json)
add_subdirectory(reference)

add_test(
  NAME ThreeBusClassical_validation
  COMMAND DynamicSimulation ThreeBusClassical.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

gridkit_example_add_file(${CMAKE_SOURCE_DIR}/cases/PhasorDynamics/Toy/ThreeBusClassical.state.json)
