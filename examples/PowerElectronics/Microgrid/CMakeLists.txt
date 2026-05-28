
add_executable(microgrid Microgrid.cpp)
target_link_libraries(
  microgrid
  GridKit::power_elec_disgen
  GridKit::power_elec_microline
  GridKit::power_elec_microload
  GridKit::solvers_dyn
  GridKit::power_elec_microbusdq)

add_test(NAME Microgrid COMMAND $<TARGET_FILE:microgrid>)
install(TARGETS microgrid RUNTIME DESTINATION bin)
