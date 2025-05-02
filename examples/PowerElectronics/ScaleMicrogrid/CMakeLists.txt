



add_executable(scalemicrogrid ScaleMicrogrid.cpp)
add_executable(scalemicrogridarbitrary ScaleMicrogridArbitrary.cpp)
target_link_libraries(scalemicrogrid GRIDKIT::power_elec_disgen
                                     GRIDKIT::power_elec_microline
                                     GRIDKIT::power_elec_microload
                                     GRIDKIT::solvers_dyn
                                     GRIDKIT::power_elec_microbusdq)

target_link_libraries(scalemicrogridarbitrary GRIDKIT::power_elec_disgen
                                        GRIDKIT::power_elec_microline
                                        GRIDKIT::power_elec_microload
                                        GRIDKIT::solvers_dyn
                                        GRIDKIT::power_elec_microbusdq)
								
add_test(NAME ScaleMicrogrid    COMMAND $<TARGET_FILE:scalemicrogrid>)
install(TARGETS scalemicrogrid RUNTIME DESTINATION bin)
