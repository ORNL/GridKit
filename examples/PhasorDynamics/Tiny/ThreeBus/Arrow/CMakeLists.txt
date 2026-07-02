gridkit_example_add_file(ThreeBusArrow.case.json)
gridkit_example_add_file(ThreeBusArrow.solver.json)

add_test(
  NAME ThreeBusArrow_pdsim
  COMMAND DynamicSimulation ThreeBusArrow.solver.json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# The Arrow IPC file format starts with the magic bytes "ARROW1"
add_test(
  NAME ThreeBusArrow_magic
  COMMAND bash -c "test -s ThreeBusArrow.arrow && head -c 6 ThreeBusArrow.arrow | grep -q ARROW1"
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(ThreeBusArrow_pdsim PROPERTIES FIXTURES_SETUP arrow_monitor_output)
set_tests_properties(ThreeBusArrow_magic PROPERTIES FIXTURES_REQUIRED arrow_monitor_output)
