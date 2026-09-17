add_executable(test_pd_tengen_classical TenGenClassical.cpp)
target_link_libraries(
  test_pd_tengen_classical
  PRIVATE GridKit::phasor_dynamics_systemmodel GridKit::solvers_dyn GridKit::testing)

configure_file(reference/TenGenClassical.ref.csv reference/TenGenClassical.ref.csv COPYONLY)

add_test(NAME TenGenClassical COMMAND test_pd_tengen_classical)

install(TARGETS test_pd_tengen_classical)
