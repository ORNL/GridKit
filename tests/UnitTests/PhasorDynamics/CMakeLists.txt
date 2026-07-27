# [[
#  Author(s):
#    - Slaven Peles <peless@ornl.gov>
#]]

add_executable(test_phasor_bus runBusTests.cpp)
target_link_libraries(
  test_phasor_bus GridKit::phasor_dynamics_bus GridKit::testing)

add_executable(test_phasor_bus_fault runBusFaultTests.cpp)
target_link_libraries(
  test_phasor_bus_fault
  GridKit::definitions
  GridKit::phasor_dynamics_bus_fault
  GridKit::phasor_dynamics_bus_fault_dependency_tracking
  GridKit::phasor_dynamics_bus
  GridKit::phasor_dynamics_bus_dependency_tracking
  GridKit::testing)

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

add_executable(test_phasor_loadz runLoadZTests.cpp)
target_link_libraries(
  test_phasor_loadz
  GridKit::definitions
  GridKit::phasor_dynamics_loadz
  GridKit::phasor_dynamics_loadz_dependency_tracking
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

add_executable(test_phasor_governor_gastpti runGovernorGastPtiTests.cpp)
target_link_libraries(
  test_phasor_governor_gastpti
  GridKit::definitions
  GridKit::phasor_dynamics_governor_gastpti
  GridKit::phasor_dynamics_governor_gastpti_dependency_tracking
  GridKit::testing)

add_executable(test_phasor_governor_hygov runGovernorHygovTests.cpp)
target_link_libraries(
  test_phasor_governor_hygov
  GridKit::definitions
  GridKit::phasor_dynamics_governor_hygov
  GridKit::phasor_dynamics_governor_hygov_dependency_tracking
  GridKit::testing)

add_executable(test_phasor_exciter_ieeet1 runExciterIeeet1Tests.cpp)
target_link_libraries(
  test_phasor_exciter_ieeet1
  GridKit::definitions
  GridKit::phasor_dynamics_components
  GridKit::phasor_dynamics_components_dependency_tracking
  GridKit::testing)

add_executable(test_phasor_exciter_esdc1a runExciterEsdc1aTests.cpp)
target_link_libraries(
  test_phasor_exciter_esdc1a
  GridKit::definitions
  GridKit::phasor_dynamics_exciter_esdc1a
  GridKit::phasor_dynamics_exciter_esdc1a_dependency_tracking
  GridKit::phasor_dynamics_bus
  GridKit::phasor_dynamics_bus_dependency_tracking
  GridKit::testing)

add_executable(test_phasor_exciter_sexspti runExciterSexsPtiTests.cpp)
target_link_libraries(
  test_phasor_exciter_sexspti
  GridKit::definitions
  GridKit::phasor_dynamics_systemmodel
  GridKit::phasor_dynamics_systemmodel_dependency_tracking
  GridKit::testing)

add_executable(test_phasor_converter_regca runConverterRegcaTests.cpp)
target_link_libraries(
  test_phasor_converter_regca
  GridKit::definitions
  GridKit::phasor_dynamics_converter_regca
  GridKit::phasor_dynamics_converter_regca_dependency_tracking
  GridKit::phasor_dynamics_bus
  GridKit::phasor_dynamics_bus_dependency_tracking
  GridKit::testing)

add_executable(test_phasor_converter_reecb runConverterReecbTests.cpp)
target_link_libraries(
  test_phasor_converter_reecb
  GridKit::definitions
  GridKit::phasor_dynamics_converter_reecb
  GridKit::phasor_dynamics_converter_reecb_dependency_tracking
  GridKit::phasor_dynamics_bus
  GridKit::phasor_dynamics_bus_dependency_tracking
  GridKit::testing)

add_executable(test_phasor_controller_repca runControllerRepcaTests.cpp)
target_link_libraries(
  test_phasor_controller_repca
  GridKit::definitions
  GridKit::phasor_dynamics_controller_repca
  GridKit::phasor_dynamics_controller_repca_dependency_tracking
  GridKit::phasor_dynamics_bus
  GridKit::phasor_dynamics_bus_dependency_tracking
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
  GridKit::phasor_dynamics_systemmodel
  GridKit::phasor_dynamics_systemmodel_dependency_tracking
  GridKit::testing)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/ThreeBusBasicBad.json
               ${CMAKE_CURRENT_BINARY_DIR}/ThreeBusBasicBad.json COPYONLY)

add_executable(test_phasor_system_single_component runSystemSingleComponentTests.cpp)
target_link_libraries(
  test_phasor_system_single_component
  GridKit::phasor_dynamics_systemmodel
  GridKit::phasor_dynamics_systemmodel_dependency_tracking
  GridKit::testing)

add_executable(test_phasor_component_connection runComponentConnectionTests.cpp)
target_link_libraries(
  test_phasor_component_connection
  GridKit::definitions
  GridKit::phasor_dynamics_systemmodel
  GridKit::testing)

add_test(NAME PhasorDynamicsBusTest COMMAND test_phasor_bus)
add_test(NAME PhasorDynamicsBusFaultTest COMMAND test_phasor_bus_fault)
add_test(NAME PhasorDynamicsBusToSignalAdapterTest COMMAND test_phasor_bustosignaladapter)
add_test(NAME PhasorDynamicsBranchTest COMMAND test_phasor_branch)
add_test(NAME PhasorDynamicsGenrouTest COMMAND test_phasor_genrou)
add_test(NAME PhasorDynamicsGovernorTgov1Test COMMAND test_phasor_governor_tgov1)
add_test(NAME PhasorDynamicsGovernorGastPtiTest COMMAND test_phasor_governor_gastpti)
add_test(NAME PhasorDynamicsGovernorHygovTest COMMAND test_phasor_governor_hygov)
add_test(NAME PhasorDynamicsExciterIeeet1Test COMMAND test_phasor_exciter_ieeet1)
add_test(NAME PhasorDynamicsExciterEsdc1aTest COMMAND test_phasor_exciter_esdc1a)
add_test(NAME PhasorDynamicsGensalTest COMMAND test_phasor_gensal)
add_test(NAME PhasorDynamicsExciterSexsPtiTest COMMAND test_phasor_exciter_sexspti)
add_test(NAME PhasorDynamicsConverterRegcaTest COMMAND test_phasor_converter_regca)
add_test(NAME PhasorDynamicsConverterReecbTest COMMAND test_phasor_converter_reecb)
add_test(NAME PhasorDynamicsControllerRepcaTest COMMAND test_phasor_controller_repca)
add_test(NAME PhasorDynamicsStabilizerIeeestTest COMMAND test_phasor_stabilizer_ieeest)
add_test(NAME PhasorDynamicsGenClassicalTest COMMAND test_phasor_gen_classical)
add_test(NAME PhasorDynamicsLoadZTest COMMAND test_phasor_loadz)
add_test(NAME PhasorDynamicsLoadZIPTest COMMAND test_phasor_loadzip)
add_test(
  NAME PhasorDynamicsSystemTest
  COMMAND test_phasor_system
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
add_test(NAME PhasorDynamicsSystemSingleComponentTest COMMAND test_phasor_system_single_component)
add_test(NAME PhasorDynamicsComponentConnectionTest COMMAND test_phasor_component_connection)

install(
  TARGETS test_phasor_bus
          test_phasor_bus_fault
          test_phasor_bustosignaladapter
          test_phasor_branch
          test_phasor_loadz
          test_phasor_loadzip
          test_phasor_genrou
          test_phasor_governor_tgov1
          test_phasor_governor_gastpti
          test_phasor_governor_hygov
          test_phasor_exciter_ieeet1
          test_phasor_exciter_esdc1a
          test_phasor_gensal
          test_phasor_exciter_sexspti
          test_phasor_converter_regca
          test_phasor_converter_reecb
          test_phasor_controller_repca
          test_phasor_stabilizer_ieeest
          test_phasor_gen_classical
          test_phasor_system
          test_phasor_system_single_component
          test_phasor_component_connection
  RUNTIME DESTINATION bin)
