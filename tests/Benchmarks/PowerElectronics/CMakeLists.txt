
add_executable(bench-rlcircuit RLCircuit.cpp)
target_link_libraries(bench-rlcircuit GridKit::power_elec_capacitor
                                GridKit::power_elec_inductor
                                GridKit::power_elec_resistor
                                GridKit::power_elec_voltagesource
                                GridKit::solvers_dyn)
                                
install(TARGETS bench-rlcircuit RUNTIME DESTINATION bin)

add_executable(bench-microgrid Microgrid.cpp)
target_link_libraries(bench-microgrid GridKit::power_elec_disgen
                                GridKit::power_elec_microline
                                GridKit::power_elec_microload
                                GridKit::solvers_dyn
                                GridKit::power_elec_microbusdq)
                                
install(TARGETS bench-microgrid RUNTIME DESTINATION bin)
