



add_executable(microgrid Microgrid.cpp)
target_link_libraries(microgrid GRIDKIT::powerelec_disgen
                                GRIDKIT::powerelec_microline
                                GRIDKIT::powerelec_microload
                                GRIDKIT::solvers_dyn
                                GRIDKIT::powerelec_microbusdq)
                                
add_test(NAME Microgrid    COMMAND $<TARGET_FILE:microgrid>)
install(TARGETS microgrid RUNTIME DESTINATION bin)
