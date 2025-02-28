add_executable(rlcircuit RLCircuit.cpp)
target_link_libraries(
  rlcircuit GRIDKIT::power_elec_capacitor GRIDKIT::power_elec_inductor
  GRIDKIT::power_elec_resistor GRIDKIT::power_elec_voltagesource
  GRIDKIT::solvers_dyn)

add_test(NAME RLCircuit COMMAND $<TARGET_FILE:rlcircuit>)
install(TARGETS rlcircuit RUNTIME DESTINATION bin)
