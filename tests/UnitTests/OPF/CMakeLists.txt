add_executable(test_opf_system_data runSystemDataTests.cpp)

target_link_libraries(
  test_opf_system_data
  PRIVATE GridKit::model_opf GridKit::testing)

add_executable(test_opf_bus_generator runBusGeneratorTests.cpp)
target_link_libraries(
  test_opf_bus_generator
  PRIVATE GridKit::model_opf GridKit::testing)

add_executable(test_opf_load_shunt_branch runLoadShuntBranchTests.cpp)
target_link_libraries(
  test_opf_load_shunt_branch
  PRIVATE GridKit::model_opf GridKit::testing)

add_executable(test_opf_system runSystemTests.cpp)
target_link_libraries(
  test_opf_system
  PRIVATE GridKit::model_opf GridKit::testing)

add_executable(test_opf_state_data_writer runStateDataWriterTests.cpp)
target_link_libraries(
  test_opf_state_data_writer
  PRIVATE GridKit::model_state_data GridKit::testing)

foreach(
  opf_case
  TwoBusBasic.opf.json
  TwoBusGensal.opf.json
  ThreeBusBasic.opf.json
  texas.opf.json)
  configure_file(${PROJECT_SOURCE_DIR}/cases/OPF/${opf_case}
                 ${CMAKE_CURRENT_BINARY_DIR}/${opf_case} COPYONLY)
endforeach()

add_test(NAME OPFSystemDataTest COMMAND test_opf_system_data)
add_test(NAME OPFBusGeneratorTest COMMAND test_opf_bus_generator)
add_test(NAME OPFLoadShuntBranchTest COMMAND test_opf_load_shunt_branch)
add_test(NAME OPFSystemTest COMMAND test_opf_system)
add_test(NAME OPFStateDataWriterTest COMMAND test_opf_state_data_writer)
set_tests_properties(OPFSystemDataTest PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

install(TARGETS test_opf_system_data)
