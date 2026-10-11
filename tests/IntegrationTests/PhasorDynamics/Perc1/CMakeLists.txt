add_executable(test_pd_perc1 Perc1IntegrationTests.cpp)
target_link_libraries(
  test_pd_perc1
  PRIVATE GridKit::phasor_dynamics_systemmodel GridKit::solvers_dyn GridKit::testing)

add_test(NAME PhasorDynamicsPerc1IntegrationTest COMMAND test_pd_perc1)
install(TARGETS test_pd_perc1)
