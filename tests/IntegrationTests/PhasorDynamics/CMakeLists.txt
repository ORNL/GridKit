add_executable(test_pd_integration runPDIntegrationTests.cpp)
target_link_libraries(
  test_pd_integration
  PRIVATE GridKit::phasor_dynamics_systemmodel GridKit::solvers_dyn GridKit::testing)

add_test(
  NAME PhasorDynamicsIntegrationTest
  COMMAND test_pd_integration
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/examples/PhasorDynamics/Tiny)

add_executable(test_pd_component_connection runComponentConnectionTests.cpp)
target_link_libraries(
  test_pd_component_connection
  PRIVATE GridKit::phasor_dynamics_systemmodel GridKit::testing)

add_test(NAME PhasorDynamicsComponentConnectionTest COMMAND test_pd_component_connection)

install(TARGETS test_pd_integration test_pd_component_connection)
