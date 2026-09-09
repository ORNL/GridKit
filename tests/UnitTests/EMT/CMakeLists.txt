# [[
#  Author(s):
#    - Luke Lowery <lukel@tamu.edu>
#]]

if(GRIDKIT_ENABLE_ENZYME AND GRIDKIT_ENABLE_SUNDIALS_SPARSE)
  add_executable(test_emt_events runEventTests.cpp)
  target_include_directories(
    test_emt_events PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/nlohmann-json/include
                            ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include)
  target_link_libraries(
    test_emt_events
    GridKit::emt_systemmodel
    GridKit::solvers_dyn
    GridKit::testing)
  add_test(NAME EMTEventsTest COMMAND test_emt_events)

  add_executable(test_emt_dae_analysis runDaeAnalysisTests.cpp)
  target_include_directories(test_emt_dae_analysis
                             PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/nlohmann-json/include)
  target_link_libraries(
    test_emt_dae_analysis
    GridKit::emt_systemmodel
    GridKit::solvers_dyn
    GridKit::testing)
  add_test(NAME EMTDaeAnalysisTest COMMAND test_emt_dae_analysis)
endif()

add_executable(test_emt_enzyme_spike runEnzymeSpikeTests.cpp)
target_link_libraries(
  test_emt_enzyme_spike
  GridKit::definitions
  GridKit::emt_bus
  GridKit::emt_bus_dependency_tracking
  GridKit::emt_voltage_source
  GridKit::emt_voltage_source_dependency_tracking
  GridKit::emt_loadz
  GridKit::emt_loadz_dependency_tracking
  GridKit::testing)
add_test(NAME EMTEnzymeSpikeTest COMMAND test_emt_enzyme_spike)

add_executable(test_emt_dependent_voltage_source runDependentVoltageSourceTests.cpp)
target_link_libraries(
  test_emt_dependent_voltage_source
  GridKit::definitions
  GridKit::emt_systemmodel
  GridKit::emt_dependent_voltage_source_dependency_tracking
  GridKit::emt_bus_dependency_tracking
  GridKit::testing)
add_test(NAME EMTDependentVoltageSourceTest COMMAND test_emt_dependent_voltage_source)

add_executable(test_emt_admittance_components runAdmittanceComponentsTests.cpp)
target_link_libraries(
  test_emt_admittance_components GridKit::emt_systemmodel GridKit::testing)
add_test(NAME EMTAdmittanceComponentsTest COMMAND test_emt_admittance_components)
install(TARGETS test_emt_admittance_components RUNTIME DESTINATION bin)

add_executable(test_emt_line_lumped runLineLumpedTests.cpp)
target_link_libraries(
  test_emt_line_lumped
  GridKit::definitions
  GridKit::emt_line_lumped
  GridKit::emt_line_lumped_dependency_tracking
  GridKit::emt_bus
  GridKit::emt_bus_dependency_tracking
  GridKit::testing)
add_test(NAME EMTLineLumpedTest COMMAND test_emt_line_lumped)

if(GRIDKIT_ENABLE_ENZYME)
  add_executable(test_emt_line_distributed runLineDistributedTests.cpp)
  target_link_libraries(test_emt_line_distributed GridKit::emt_systemmodel GridKit::testing)
  add_test(NAME EMTLineDistributedTest COMMAND test_emt_line_distributed)
endif()

add_executable(test_emt_machine runMachineTests.cpp)
target_link_libraries(
  test_emt_machine
  GridKit::definitions
  GridKit::emt_machine
  GridKit::emt_machine_dependency_tracking
  GridKit::emt_bus
  GridKit::emt_bus_dependency_tracking
  GridKit::testing)
add_test(NAME EMTMachineTest COMMAND test_emt_machine)

add_executable(test_emt_transformer runTransformerTests.cpp)
target_link_libraries(
  test_emt_transformer
  GridKit::definitions
  GridKit::emt_transformer
  GridKit::emt_transformer_dependency_tracking
  GridKit::emt_bus
  GridKit::emt_bus_dependency_tracking
  GridKit::testing)
add_test(NAME EMTTransformerTest COMMAND test_emt_transformer)

add_executable(test_emt_pll runPllTests.cpp)
target_link_libraries(
  test_emt_pll
  GridKit::definitions
  GridKit::emt_pll
  GridKit::emt_pll_dependency_tracking
  GridKit::emt_bus
  GridKit::emt_bus_dependency_tracking
  GridKit::testing)
add_test(NAME EMTPllTest COMMAND test_emt_pll)

add_executable(test_emt_filter runFilterTests.cpp)
target_include_directories(
  test_emt_filter
  PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/nlohmann-json/include
          ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include)
target_link_libraries(
  test_emt_filter
  GridKit::emt_systemmodel
  GridKit::emt_filter_dependency_tracking
  GridKit::testing)
add_test(NAME EMTFilterTest COMMAND test_emt_filter)

add_executable(test_emt_repca runRepcaTests.cpp)
target_link_libraries(
  test_emt_repca
  GridKit::definitions
  GridKit::emt_repca
  GridKit::phasor_dynamics_controller_repca
  GridKit::testing)
add_test(NAME EMTRepcaTest COMMAND test_emt_repca)

add_executable(test_emt_reecb runReecbTests.cpp)
target_link_libraries(
  test_emt_reecb
  GridKit::definitions
  GridKit::emt_reecb
  GridKit::phasor_dynamics_controller_reecb
  GridKit::testing)
add_test(NAME EMTReecbTest COMMAND test_emt_reecb)

add_executable(test_emt_inner_current_control runInnerCurrentControlTests.cpp)
target_link_libraries(
  test_emt_inner_current_control
  GridKit::definitions
  GridKit::emt_inner_current_control
  GridKit::emt_inner_current_control_dependency_tracking
  GridKit::testing)
add_test(NAME EMTInnerCurrentControlTest COMMAND test_emt_inner_current_control)

add_executable(test_emt_vector_fit runVectorFitTests.cpp)
target_link_libraries(
  test_emt_vector_fit
  GridKit::definitions
  GridKit::emt_vector_fit
  GridKit::emt_vector_fit_dependency_tracking
  GridKit::testing)
add_test(NAME EMTVectorFitTest COMMAND test_emt_vector_fit)

add_executable(test_emt_system runSystemTests.cpp)
target_include_directories(
  test_emt_system
  PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/nlohmann-json/include)
target_link_libraries(
  test_emt_system
  GridKit::definitions
  GridKit::emt_systemmodel
  GridKit::testing)
if(GRIDKIT_ENABLE_SUNDIALS)
  target_link_libraries(test_emt_system GridKit::solvers_dyn)
  target_compile_definitions(test_emt_system PRIVATE GRIDKIT_ENABLE_SUNDIALS)
endif()
add_test(NAME EMTSystemTest COMMAND test_emt_system)

add_executable(test_emt_current_state runCurrentStateTests.cpp)
target_include_directories(
  test_emt_current_state
  PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/nlohmann-json/include)
target_link_libraries(
  test_emt_current_state
  GridKit::emt_systemmodel
  GridKit::model_state_data
  GridKit::testing)
add_test(NAME EMTCurrentStateTest COMMAND test_emt_current_state)

add_executable(test_emt_pwm_converter runPwmConverterTests.cpp)
target_include_directories(
  test_emt_pwm_converter
  PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/nlohmann-json/include
          ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include)
target_link_libraries(
  test_emt_pwm_converter
  GridKit::emt_systemmodel
  GridKit::emt_pwm_dependency_tracking
  GridKit::emt_converter_dependency_tracking
  GridKit::testing)
if(GRIDKIT_ENABLE_SUNDIALS)
  target_link_libraries(test_emt_pwm_converter GridKit::solvers_dyn)
  target_compile_definitions(test_emt_pwm_converter PRIVATE GRIDKIT_ENABLE_SUNDIALS)
endif()
add_test(NAME EMTPwmConverterTest COMMAND test_emt_pwm_converter)

install(
  TARGETS test_emt_pwm_converter
          test_emt_dependent_voltage_source
          test_emt_enzyme_spike
          test_emt_line_lumped
          test_emt_machine
          test_emt_transformer
          test_emt_system
          test_emt_vector_fit
  RUNTIME DESTINATION bin)

if(GRIDKIT_ENABLE_SUNDIALS)
  add_executable(test_emt_study_data runStudyDataTests.cpp)
  target_include_directories(
    test_emt_study_data PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/nlohmann-json/include
                                ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include)
  target_link_libraries(
    test_emt_study_data
    GridKit::emt_systemmodel
    GridKit::solvers_dyn
    GridKit::testing)
  add_test(NAME EMTStudyDataTest COMMAND test_emt_study_data)
endif()

add_executable(test_emt_ieeet1 runIeeet1Tests.cpp)
target_link_libraries(
  test_emt_ieeet1
  GridKit::definitions
  GridKit::emt_ieeet1
  GridKit::emt_ieeet1_dependency_tracking
  GridKit::testing)
add_test(NAME EMTIeeet1Test COMMAND test_emt_ieeet1)
install(TARGETS test_emt_ieeet1 RUNTIME DESTINATION bin)

add_executable(test_emt_rational runRationalTests.cpp)
target_include_directories(
  test_emt_rational
  PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/nlohmann-json/include)
target_link_libraries(
  test_emt_rational
  GridKit::emt_vector_fit
  GridKit::emt_state_space
  GridKit::testing)
add_test(NAME EMTRationalTest COMMAND test_emt_rational)

add_executable(test_emt_sexs_pti runSexsPtiTests.cpp)
target_link_libraries(
  test_emt_sexs_pti
  GridKit::emt_sexs_pti
  GridKit::emt_sexs_pti_dependency_tracking
  GridKit::phasor_dynamics_exciter_sexspti
  GridKit::phasor_dynamics_bus
  GridKit::testing)
add_test(NAME EMTSexsPtiTest COMMAND test_emt_sexs_pti)
install(TARGETS test_emt_sexs_pti RUNTIME DESTINATION bin)

add_executable(test_emt_gastpti runGastPtiTests.cpp)
target_link_libraries(
  test_emt_gastpti
  GridKit::emt_gastpti
  GridKit::emt_gastpti_dependency_tracking
  GridKit::phasor_dynamics_governor_gastpti
  GridKit::testing)
add_test(NAME EMTGastPtiTest COMMAND test_emt_gastpti)
install(TARGETS test_emt_gastpti RUNTIME DESTINATION bin)

add_executable(test_emt_ieeest runIeeestTests.cpp)
target_link_libraries(
  test_emt_ieeest
  GridKit::emt_ieeest
  GridKit::emt_ieeest_dependency_tracking
  GridKit::phasor_dynamics_stabilizer_ieeest
  GridKit::testing)
add_test(NAME EMTIeeestTest COMMAND test_emt_ieeest)
install(TARGETS test_emt_ieeest RUNTIME DESTINATION bin)

add_executable(test_emt_model_data runModelDataTests.cpp)
target_include_directories(
  test_emt_model_data
  PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/nlohmann-json/include
          ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include)
target_link_libraries(
  test_emt_model_data GridKit::emt_systemmodel GridKit::testing)
add_test(NAME EMTModelDataTest COMMAND test_emt_model_data)

add_executable(test_emt_initialization runInitializationTests.cpp)
target_include_directories(
  test_emt_initialization
  PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/nlohmann-json/include
          ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include)
target_link_libraries(
  test_emt_initialization GridKit::emt_systemmodel GridKit::testing)
add_test(NAME EMTInitializationTest COMMAND test_emt_initialization)

add_executable(test_emt_regfma runRegfmaTests.cpp)
target_link_libraries(
  test_emt_regfma
  GridKit::emt_systemmodel
  GridKit::emt_regfma_dependency_tracking
  GridKit::testing)
if(GRIDKIT_ENABLE_ENZYME AND GRIDKIT_ENABLE_SUNDIALS_SPARSE)
  target_link_libraries(test_emt_regfma GridKit::solvers_dyn)
  target_compile_definitions(test_emt_regfma PRIVATE GRIDKIT_ENABLE_SUNDIALS)
endif()
add_test(NAME EMTRegfmaTest COMMAND test_emt_regfma)
install(TARGETS test_emt_regfma RUNTIME DESTINATION bin)

add_executable(test_variable_monitor runVariableMonitorTests.cpp)
target_include_directories(
  test_variable_monitor
  PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/nlohmann-json/include
          ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include)
target_link_libraries(
  test_variable_monitor GridKit::definitions GridKit::testing)
add_test(NAME VariableMonitorTest COMMAND test_variable_monitor)
