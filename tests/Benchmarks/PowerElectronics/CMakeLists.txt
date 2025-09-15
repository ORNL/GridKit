
add_executable(bench-rlcircuit RLCircuit.cpp)
target_link_libraries(bench-rlcircuit GridKit::power_elec_capacitor
                                GridKit::power_elec_inductor
                                GridKit::power_elec_resistor
                                GridKit::power_elec_voltagesource
                                GridKit::solvers_dyn)
                                
install(TARGETS bench-rlcircuit RUNTIME DESTINATION bin)
