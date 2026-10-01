add_executable(microgrid Microgrid.cpp)

target_link_libraries(
  microgrid GridKit::solvers_dyn power_elec_microgrid_network)

add_test(NAME Microgrid COMMAND $<TARGET_FILE:microgrid>)
install(TARGETS microgrid RUNTIME DESTINATION bin)
