add_executable(test_pd_system_jacobian runSystemJacobianTests.cpp)
target_link_libraries(
  test_pd_system_jacobian
  PRIVATE GridKit::phasor_dynamics_systemmodel
          GridKit::phasor_dynamics_systemmodel_dependency_tracking
          GridKit::utilities_cli_args
          GridKit::testing)

# All cases in cases/PhasorDynamics, one sub-test each. One evaluation because
# this checks the Jacobians, meaningful timings need the driver's default.
add_test(
  NAME PhasorDynamicsIntegrationSystemJacobian
  COMMAND test_pd_system_jacobian --cases ${CMAKE_SOURCE_DIR}/cases/PhasorDynamics --evaluations 1)

install(TARGETS test_pd_system_jacobian)
