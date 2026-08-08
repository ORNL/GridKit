# The case file stays in cases/PhasorDynamics, so the solver file is run from
# the source tree; system_model_file is relative to the solver file's location.
add_test(
  NAME socar_validation
  COMMAND DynamicSimulation ${CMAKE_CURRENT_SOURCE_DIR}/socar.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
