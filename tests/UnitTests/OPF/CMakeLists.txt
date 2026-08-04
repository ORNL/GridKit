add_executable(test_opf_state_data runStateDataTests.cpp)
target_link_libraries(
  test_opf_state_data
  PRIVATE GridKit::model_state_data GridKit::testing)
add_test(NAME OPFStateDataTest COMMAND test_opf_state_data)

if(TARGET model_opf)
  add_executable(test_opf_components runComponentTests.cpp)
  target_link_libraries(test_opf_components PRIVATE GridKit::model_opf_components GridKit::testing)
  add_test(NAME OPFComponentTest COMMAND test_opf_components)

  add_executable(test_opf_system_data_parser runSystemDataParserTests.cpp)
  target_link_libraries(test_opf_system_data_parser PRIVATE GridKit::model_opf GridKit::testing)
  target_compile_definitions(test_opf_system_data_parser
                             PRIVATE GRIDKIT_OPF_CASE_DIRECTORY="${PROJECT_SOURCE_DIR}/cases/OPF")
  add_test(NAME OPFSystemDataParserTest COMMAND test_opf_system_data_parser)

  add_executable(test_opf_system runSystemTests.cpp)
  target_link_libraries(test_opf_system PRIVATE GridKit::model_opf GridKit::testing)
  add_test(NAME OPFSystemTest COMMAND test_opf_system)
endif()
