add_executable(test_pd_integration runPDIntegrationTests.cpp)
target_link_libraries(
  test_pd_integration
  PRIVATE GridKit::phasor_dynamics_systemmodel GridKit::solvers_dyn GridKit::testing)

add_test(
  NAME PhasorDynamicsIntegrationTest
  COMMAND test_pd_integration
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/examples/PhasorDynamics/Validation/Toy)

install(TARGETS test_pd_integration)

add_executable(test_pd_signal_delay runSignalDelayIntegrationTests.cpp)
target_link_libraries(
  test_pd_signal_delay
  PRIVATE GridKit::phasor_dynamics_systemmodel GridKit::solvers_dyn GridKit::testing)
add_test(NAME PhasorDynamicsSignalDelayIntegrationTest COMMAND test_pd_signal_delay)
set_tests_properties(PhasorDynamicsSignalDelayIntegrationTest PROPERTIES SKIP_RETURN_CODE 77)
install(TARGETS test_pd_signal_delay)

add_subdirectory(TenGen)
