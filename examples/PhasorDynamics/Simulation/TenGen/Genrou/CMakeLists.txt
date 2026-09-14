add_executable(TenGenGenrou TenGenGenrou.cpp)
target_link_libraries(
  TenGenGenrou GridKit::phasor_dynamics_systemmodel GridKit::solvers_dyn)
install(TARGETS TenGenGenrou RUNTIME DESTINATION bin)

add_test(NAME TenGenGenrou COMMAND $<TARGET_FILE:TenGenGenrou>)
