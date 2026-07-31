add_executable(test_opf_system_data runSystemDataTests.cpp)

target_link_libraries(
  test_opf_system_data
  PRIVATE GridKit::model_opf GridKit::testing)

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
set_tests_properties(OPFSystemDataTest PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

install(TARGETS test_opf_system_data)
