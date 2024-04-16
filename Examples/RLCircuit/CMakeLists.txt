



add_executable(rlcircuit RLCircuit.cpp)
target_link_libraries(rlcircuit GRIDKIT::powerelec_capacitor
                                GRIDKIT::powerelec_inductor
                                GRIDKIT::powerelec_resistor
                                GRIDKIT::powerelec_voltagesource
                                GRIDKIT::solvers_dyn)
                                
add_test(NAME RLCircuit    COMMAND $<TARGET_FILE:rlcircuit>)
install(TARGETS rlcircuit RUNTIME DESTINATION bin)
