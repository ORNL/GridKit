gridkit_example_add_file(${CMAKE_SOURCE_DIR}/cases/PhasorDynamics/WECC240/WECC240.case.json)
gridkit_example_add_file(${CMAKE_SOURCE_DIR}/cases/PhasorDynamics/WECC240/WECC240.state.json)
gridkit_example_add_file(WECC240.solver.json)
gridkit_example_add_file(WECC240.partition.json)

add_test(
  NAME WECC240_partitioned
  COMMAND DynamicSimulation WECC240.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

if(GRIDKIT_ENABLE_SUNDIALS_PARTITIONED_DAE)
  gridkit_example_add_file(WECC240.integrator.solver.json)
  # TODO: Run the fault study once events restart from consistent algebraic
  # variables and partitions are solved directly. Until then the test runs two steps.
  add_test(
    NAME WECC240_partitioned_integrator
    COMMAND PartitionedDynamicSimulation WECC240.integrator.solver.json
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endif()
