# [[
#  Author(s):
#    - Luke Lowery <lukel@tamu.edu>
#]]

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
