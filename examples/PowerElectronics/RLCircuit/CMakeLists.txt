
add_executable(rlcircuit RLCircuit.cpp)
target_link_libraries(
  rlcircuit
  GridKit::power_elec_capacitor
  GridKit::power_elec_inductor
  GridKit::power_elec_resistor
  GridKit::power_elec_voltagesource
  GridKit::solvers_dyn)

add_test(NAME RLCircuit COMMAND $<TARGET_FILE:rlcircuit>)
install(TARGETS rlcircuit RUNTIME DESTINATION bin)
