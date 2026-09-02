# [[
#  Author(s):
#    - Luke Lowery <lukel@tamu.edu>
#]]

add_executable(test_emt_enzyme_spike runEnzymeSpikeTests.cpp)
target_link_libraries(
  test_emt_enzyme_spike
  GridKit::definitions
  GridKit::emt_bus
  GridKit::emt_bus_dependency_tracking
  GridKit::emt_voltage_source
  GridKit::emt_voltage_source_dependency_tracking
  GridKit::emt_loadz
  GridKit::emt_loadz_dependency_tracking
  GridKit::testing)
add_test(NAME EMTEnzymeSpikeTest COMMAND test_emt_enzyme_spike)

add_executable(test_emt_line_lumped runLineLumpedTests.cpp)
target_link_libraries(
  test_emt_line_lumped
  GridKit::definitions
  GridKit::emt_line_lumped
  GridKit::emt_line_lumped_dependency_tracking
  GridKit::emt_bus
  GridKit::emt_bus_dependency_tracking
  GridKit::testing)
add_test(NAME EMTLineLumpedTest COMMAND test_emt_line_lumped)

add_executable(test_emt_synchronous_machine runSynchronousMachineTests.cpp)
target_link_libraries(
  test_emt_synchronous_machine
  GridKit::definitions
  GridKit::emt_synchronous_machine
  GridKit::emt_synchronous_machine_dependency_tracking
  GridKit::emt_bus
  GridKit::emt_bus_dependency_tracking
  GridKit::testing)
add_test(NAME EMTSynchronousMachineTest COMMAND test_emt_synchronous_machine)

add_executable(test_emt_vector_fit runVectorFitTests.cpp)
target_link_libraries(
  test_emt_vector_fit
  GridKit::definitions
  GridKit::emt_vector_fit
  GridKit::emt_vector_fit_dependency_tracking
  GridKit::testing)
add_test(NAME EMTVectorFitTest COMMAND test_emt_vector_fit)

add_executable(test_emt_system runSystemTests.cpp)
target_link_libraries(
  test_emt_system
  GridKit::definitions
  GridKit::emt_systemmodel
  GridKit::testing)
if(GRIDKIT_ENABLE_SUNDIALS)
  target_link_libraries(test_emt_system GridKit::solvers_dyn)
  target_compile_definitions(test_emt_system PRIVATE GRIDKIT_ENABLE_SUNDIALS)
endif()
add_test(NAME EMTSystemTest COMMAND test_emt_system)

install(
  TARGETS test_emt_enzyme_spike
          test_emt_line_lumped
          test_emt_synchronous_machine
          test_emt_system
          test_emt_vector_fit
  RUNTIME DESTINATION bin)
