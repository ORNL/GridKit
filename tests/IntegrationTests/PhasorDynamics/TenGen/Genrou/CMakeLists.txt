add_executable(test_pd_tengen_genrou TenGenGenrou.cpp)
target_link_libraries(
  test_pd_tengen_genrou
  PRIVATE GridKit::phasor_dynamics_systemmodel GridKit::solvers_dyn GridKit::testing)

configure_file(reference/TenGenGenrou.ref.csv reference/TenGenGenrou.ref.csv COPYONLY)

add_test(NAME TenGenGenrou COMMAND test_pd_tengen_genrou)

install(TARGETS test_pd_tengen_genrou)
install(FILES reference/TenGenGenrou.ref.csv
        DESTINATION share/gridkit/tests/PhasorDynamics/TenGen/Genrou/reference)
