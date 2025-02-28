
add_executable(microgrid Microgrid.cpp)
target_link_libraries(
  microgrid GRIDKIT::power_elec_disgen GRIDKIT::power_elec_microline
  GRIDKIT::power_elec_microload GRIDKIT::solvers_dyn
  GRIDKIT::power_elec_microbusdq)

add_test(NAME Microgrid COMMAND $<TARGET_FILE:microgrid>)
install(TARGETS microgrid RUNTIME DESTINATION bin)
