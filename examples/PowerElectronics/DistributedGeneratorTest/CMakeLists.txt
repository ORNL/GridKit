add_executable(dgtest DGTest.cpp)
target_link_libraries(
  dgtest
  GridKit::power_elec_disgen
  GridKit::power_elec_microline
  GridKit::power_elec_microload)

add_test(NAME DistributedGeneratorTest COMMAND $<TARGET_FILE:dgtest>)
install(TARGETS dgtest RUNTIME DESTINATION bin)
