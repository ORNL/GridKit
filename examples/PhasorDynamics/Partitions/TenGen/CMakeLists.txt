if(GRIDKIT_ENABLE_SUNDIALS_PARTITIONED_DAE)
  gridkit_example_add_file(TenGen.case.json)
  gridkit_example_add_file(TenGen.partition.json)
  gridkit_example_add_file(TenGen.integrator.solver.json)
  add_test(
    NAME TenGen_partitioned_integrator
    COMMAND PartitionedDynamicSimulation TenGen.integrator.solver.json
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endif()
