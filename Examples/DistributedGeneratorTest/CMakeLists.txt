



add_executable(dgtest DGTest.cpp)
target_link_libraries(dgtest GRIDKIT::powerelec_disgen
                                GRIDKIT::powerelec_mircoline
                                GRIDKIT::powerelec_microload
                                GRIDKIT::solvers_dyn)
                                
add_test(NAME DistributedGeneratorTest    COMMAND $<TARGET_FILE:dgtest>)
install(TARGETS dgtest RUNTIME DESTINATION bin)
