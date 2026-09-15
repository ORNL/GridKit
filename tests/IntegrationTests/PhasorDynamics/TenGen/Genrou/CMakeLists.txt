add_executable(test_pd_tengen_genrou TenGenGenrou.cpp)
target_link_libraries(
  test_pd_tengen_genrou
  PRIVATE GridKit::phasor_dynamics_systemmodel GridKit::solvers_dyn)

add_test(NAME TenGenGenrou COMMAND test_pd_tengen_genrou)

install(TARGETS test_pd_tengen_genrou)
