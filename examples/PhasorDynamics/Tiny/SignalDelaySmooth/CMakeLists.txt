gridkit_example_add_file(SignalDelaySmooth.case.json)
gridkit_example_add_file(SignalDelaySmooth.solver.json)
gridkit_example_add_file(SignalDelaySmoothSource.csv)

add_test(
  NAME SignalDelaySmooth_pdsim
  COMMAND DynamicSimulation SignalDelaySmooth.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
