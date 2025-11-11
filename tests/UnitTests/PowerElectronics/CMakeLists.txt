add_executable(test_power_elec_system runSystemTests.cpp)
target_link_libraries(test_power_elec_system GridKit::power_elec_disgen
                                             GridKit::power_elec_microline
                                             GridKit::power_elec_microload
                                             GridKit::power_elec_microbusdq)

add_test(NAME PowerElectronicsSystemTest COMMAND test_power_elec_system)

install(TARGETS test_power_elec_system
        RUNTIME DESTINATION bin)
