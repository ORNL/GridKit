



add_executable(dgtest DGTest.cpp)
target_link_libraries(dgtest GRIDKIT::powerelec_disgen)
								
add_test(NAME DiscreteGeneratorTest    COMMAND $<TARGET_FILE:dgtest>)
install(TARGETS dgtest RUNTIME DESTINATION bin)
