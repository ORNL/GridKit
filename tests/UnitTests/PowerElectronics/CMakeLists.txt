add_executable(test_power_electronics_node runCircuitNodeTests.cpp)
target_link_libraries(
  test_power_electronics_node
  PRIVATE GridKit::power_electronics_circuit_node GridKit::testing)

add_executable(test_subsystem_model_with_hires runSubsystemModelWithHiresTest.cpp)
target_link_libraries(
  test_subsystem_model_with_hires
  PRIVATE GridKit::power_elec_partition_interfaces GridKit::testing GridKit::sparse_matrix)

add_executable(test_power_electronics_component_clone runComponentCloneTests.cpp)
target_link_libraries(
  test_power_electronics_component_clone
  PRIVATE GridKit::power_elec_disgen
          GridKit::power_elec_microline
          GridKit::power_elec_microload
          GridKit::power_elec_microbusdq
          GridKit::testing)

add_test(NAME PowerElectronicsNodeTest COMMAND $<TARGET_FILE:test_power_electronics_node>)
add_test(NAME SubsystemModelWithHires COMMAND $<TARGET_FILE:test_subsystem_model_with_hires>)
add_test(NAME PowerElectronicsComponentCloneTest
         COMMAND $<TARGET_FILE:test_power_electronics_component_clone>)

install(TARGETS test_power_electronics_node RUNTIME DESTINATION bin)
install(TARGETS test_subsystem_model_with_hires RUNTIME DESTINATION bin)
install(TARGETS test_power_electronics_component_clone RUNTIME DESTINATION bin)
