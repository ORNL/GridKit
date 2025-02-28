# [[ Author(s): - Slaven Peles <peless@ornl.gov> ]]

add_executable(test_sparsity_pattern runSparsityPatternTests.cpp)
target_link_libraries(test_sparsity_pattern GRIDKIT::phasor_dynamics_bus)

add_test(NAME SparsityPatternTest COMMAND $<TARGET_FILE:test_sparsity_pattern>)

install(TARGETS test_sparsity_pattern RUNTIME DESTINATION bin)
