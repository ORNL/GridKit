add_executable(test_pd_integration runPDIntegrationTests.cpp)
target_link_libraries(
  test_pd_integration
  PRIVATE GridKit::phasor_dynamics_systemmodel GridKit::solvers_dyn GridKit::testing)

add_test(
  NAME PhasorDynamicsIntegrationTest
  COMMAND test_pd_integration
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/examples/PhasorDynamics/Tiny)

install(TARGETS test_pd_integration)
