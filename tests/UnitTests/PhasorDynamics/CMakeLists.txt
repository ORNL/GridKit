# [[
#  Author(s):
#    - Slaven Peles <peless@ornl.gov>
#]]

add_executable(test_phasor_bus runBusTests.cpp)
target_link_libraries(
  test_phasor_bus GridKit::phasor_dynamics_bus GridKit::testing)

add_executable(test_phasor_bustosignaladapter runBusToSignalAdapterTests.cpp)
target_link_libraries(
  test_phasor_bustosignaladapter
  GridKit::phasor_dynamics_bus
  GridKit::phasor_dynamics_bus_dependency_tracking
  GridKit::phasor_dynamics_bustosignaladapter
  GridKit::phasor_dynamics_bustosignaladapter_dependency_tracking
  GridKit::testing)

add_executable(test_phasor_branch runBranchTests.cpp)
target_link_libraries(
  test_phasor_branch
  GridKit::phasor_dynamics_branch
  GridKit::phasor_dynamics_branch_dependency_tracking
  GridKit::phasor_dynamics_bus
  GridKit::phasor_dynamics_bus_dependency_tracking
  GridKit::testing)

add_executable(test_phasor_load runLoadTests.cpp)
target_link_libraries(
  test_phasor_load
  GridKit::definitions
  GridKit::phasor_dynamics_load
  GridKit::phasor_dynamics_load_dependency_tracking
  GridKit::phasor_dynamics_bus
  GridKit::phasor_dynamics_bus_dependency_tracking
  GridKit::testing)

add_executable(test_phasor_loadzip runLoadZIPTests.cpp)
target_link_libraries(
  test_phasor_loadzip
  GridKit::definitions
  GridKit::phasor_dynamics_loadzip
  GridKit::phasor_dynamics_loadzip_dependency_tracking
  GridKit::phasor_dynamics_bus
  GridKit::phasor_dynamics_bus_dependency_tracking
  GridKit::testing)

add_executable(test_phasor_genrou runGenrouTests.cpp)
target_link_libraries(
  test_phasor_genrou
  GridKit::definitions
  GridKit::phasor_dynamics_genrou
  GridKit::phasor_dynamics_genrou_dependency_tracking
  GridKit::phasor_dynamics_bus
  GridKit::phasor_dynamics_bus_dependency_tracking
  GridKit::testing)

add_executable(test_phasor_gensal runGensalTests.cpp)
target_link_libraries(
  test_phasor_gensal
  GridKit::definitions
  GridKit::phasor_dynamics_gensal
  GridKit::phasor_dynamics_gensal_dependency_tracking
  GridKit::phasor_dynamics_bus
  GridKit::phasor_dynamics_bus_dependency_tracking
  GridKit::testing)

add_executable(test_phasor_governor_tgov1 runGovernorTgov1Tests.cpp)
target_link_libraries(
  test_phasor_governor_tgov1
  GridKit::definitions
  GridKit::phasor_dynamics_genrou
  GridKit::phasor_dynamics_genrou_dependency_tracking
  GridKit::phasor_dynamics_governor_tgov1
  GridKit::phasor_dynamics_governor_tgov1_dependency_tracking
  GridKit::phasor_dynamics_signal
  GridKit::phasor_dynamics_bus
  GridKit::phasor_dynamics_bus_dependency_tracking
  GridKit::testing)

add_executable(test_phasor_exciter_ieeet1 runExciterIeeet1Tests.cpp)
target_link_libraries(
  test_phasor_exciter_ieeet1
  GridKit::definitions
  GridKit::phasor_dynamics_components
  GridKit::phasor_dynamics_components_dependency_tracking
  GridKit::testing)

add_executable(test_phasor_exciter_sexspti runExciterSexsPtiTests.cpp)
target_link_libraries(
  test_phasor_exciter_sexspti
  GridKit::definitions
  GridKit::phasor_dynamics_components
  GridKit::phasor_dynamics_components_dependency_tracking
  GridKit::testing)

add_executable(test_phasor_converter_regca runConverterRegcaTests.cpp)
target_link_libraries(
  test_phasor_converter_regca
  GridKit::definitions
  GridKit::phasor_dynamics_components
  GridKit::phasor_dynamics_components_dependency_tracking
  GridKit::testing)

add_executable(test_phasor_converter_repca runConverterRepcaTests.cpp)
target_link_libraries(
  test_phasor_converter_repca
  GridKit::definitions
  GridKit::phasor_dynamics_components
  GridKit::phasor_dynamics_components_dependency_tracking
  GridKit::testing)

add_executable(test_phasor_stabilizer_ieeest runStabilizerIeeestTests.cpp)
target_link_libraries(
  test_phasor_stabilizer_ieeest
  GridKit::definitions
  GridKit::phasor_dynamics_stabilizer_ieeest
  GridKit::phasor_dynamics_stabilizer_ieeest_dependency_tracking
  GridKit::phasor_dynamics_signal
  GridKit::testing)

add_executable(test_phasor_gen_classical runGenClassicalTests.cpp)
target_link_libraries(
  test_phasor_gen_classical
  GridKit::definitions
  GridKit::phasor_dynamics_gen_classical
  GridKit::phasor_dynamics_gen_classical_dependency_tracking
  GridKit::phasor_dynamics_bus
  GridKit::phasor_dynamics_bus_dependency_tracking
  GridKit::testing)

add_executable(test_phasor_system runSystemTests.cpp)
target_link_libraries(
  test_phasor_system
  GridKit::definitions
  GridKit::phasor_dynamics_components
  GridKit::phasor_dynamics_components_dependency_tracking
  GridKit::testing)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/ThreeBusBasicBad.json
               ${CMAKE_CURRENT_BINARY_DIR}/ThreeBusBasicBad.json COPYONLY)

add_executable(test_phasor_system_single_component runSystemSingleComponentTests.cpp)
target_link_libraries(
  test_phasor_system_single_component
  GridKit::phasor_dynamics_components
  GridKit::phasor_dynamics_components_dependency_tracking
  GridKit::testing)

add_test(NAME PhasorDynamicsBusTest COMMAND test_phasor_bus)
add_test(NAME PhasorDynamicsBusToSignalAdapterTest COMMAND test_phasor_bustosignaladapter)
add_test(NAME PhasorDynamicsBranchTest COMMAND test_phasor_branch)
add_test(NAME PhasorDynamicsGenrouTest COMMAND test_phasor_genrou)
add_test(NAME PhasorDynamicsGovernorTgov1Test COMMAND test_phasor_governor_tgov1)
add_test(NAME PhasorDynamicsExciterIeeet1Test COMMAND test_phasor_exciter_ieeet1)
add_test(NAME PhasorDynamicsGensalTest COMMAND test_phasor_gensal)
add_test(NAME PhasorDynamicsExciterSexsPtiTest COMMAND test_phasor_exciter_sexspti)
add_test(NAME PhasorDynamicsConverterRegcaTest COMMAND test_phasor_converter_regca)
add_test(NAME PhasorDynamicsConverterRepcaTest COMMAND test_phasor_converter_repca)
add_test(NAME PhasorDynamicsStabilizerIeeestTest COMMAND test_phasor_stabilizer_ieeest)
add_test(NAME PhasorDynamicsGenClassicalTest COMMAND test_phasor_gen_classical)
add_test(NAME PhasorDynamicsLoadTest COMMAND test_phasor_load)
add_test(NAME PhasorDynamicsLoadZIPTest COMMAND test_phasor_loadzip)
add_test(
  NAME PhasorDynamicsSystemTest
  COMMAND test_phasor_system
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
add_test(NAME PhasorDynamicsSystemSingleComponentTest COMMAND test_phasor_system_single_component)

install(
  TARGETS test_phasor_bus
          test_phasor_bustosignaladapter
          test_phasor_branch
          test_phasor_load
          test_phasor_loadzip
          test_phasor_genrou
          test_phasor_governor_tgov1
          test_phasor_exciter_ieeet1
          test_phasor_gensal
          test_phasor_exciter_sexspti
          test_phasor_converter_regca
          test_phasor_converter_repca
          test_phasor_stabilizer_ieeest
          test_phasor_gen_classical
          test_phasor_system
          test_phasor_system_single_component
  RUNTIME DESTINATION bin)
