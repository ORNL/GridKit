



add_executable(scalemicrogrid ScaleMicrogrid.cpp)
target_link_libraries(scalemicrogrid GRIDKIT::powerelec_disgen
                                     GRIDKIT::powerelec_microline
                                     GRIDKIT::powerelec_microload
                                     GRIDKIT::solvers_dyn
                                     GRIDKIT::powerelec_microbusdq)
								
add_test(NAME ScaleMicrogrid    COMMAND $<TARGET_FILE:scalemicrogrid>)
install(TARGETS scalemicrogrid RUNTIME DESTINATION bin)
