



add_executable(microgridfault MicrogridFault.cpp)
target_link_libraries(microgridfault GridKit::power_elec_disgen
                                     GridKit::power_elec_microline
                                     GridKit::power_elec_microload
                                     GridKit::solvers_dyn
                                     GridKit::power_elec_microbusdq)
                                
add_test(NAME MicrogridFault    COMMAND $<TARGET_FILE:microgridfault>)
install(TARGETS microgridfault RUNTIME DESTINATION bin)
