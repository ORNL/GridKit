add_executable(test_opf_steady_start runSteadyStartTests.cpp)
target_link_libraries(
  test_opf_steady_start
  PRIVATE GridKit::phasor_dynamics_systemmodel GridKit::solvers_dyn GridKit::testing)

install(TARGETS test_opf_steady_start)
