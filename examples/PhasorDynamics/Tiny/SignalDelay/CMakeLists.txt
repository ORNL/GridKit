gridkit_example_add_file(SignalDelay.case.json)
gridkit_example_add_file(SignalDelay.solver.json)
gridkit_example_add_file(SignalDelaySource.csv)

add_test(
  NAME SignalDelay_pdsim
  COMMAND DynamicSimulation SignalDelay.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
