



add_executable(dgtest DGTest.cpp)
target_link_libraries(dgtest GRIDKIT::power_elec_disgen
                                GRIDKIT::power_elec_microline
                                GRIDKIT::power_elec_microload)
                                #GRIDKIT::solvers_dyn)
                                
add_test(NAME DistributedGeneratorTest    COMMAND $<TARGET_FILE:dgtest>)
install(TARGETS dgtest RUNTIME DESTINATION bin)
