add_executable(test_power_electronics_node runCircuitNodeTests.cpp)
target_link_libraries(test_power_electronics_node PRIVATE GridKit::power_electronics_circuit_node)

add_test(NAME PowerElectronicsNodeTest COMMAND $<TARGET_FILE:test_power_electronics_node>)

install(TARGETS test_power_electronics_node
        RUNTIME DESTINATION bin)