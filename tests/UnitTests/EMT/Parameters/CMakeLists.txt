add_executable(test_emt_overhead runOverheadTests.cpp)
target_link_libraries(
  test_emt_overhead
  GridKit::emt_overhead
  GridKit::testing
  Eigen3::Eigen)

add_test(NAME EMTOverheadTest COMMAND $<TARGET_FILE:test_emt_overhead>)

add_executable(test_emt_conductor runConductorTests.cpp)
target_link_libraries(
  test_emt_conductor GridKit::emt_conductor GridKit::testing)

add_test(NAME EMTConductorTest COMMAND $<TARGET_FILE:test_emt_conductor>)

add_executable(test_emt_path runPathTests.cpp)
target_link_libraries(
  test_emt_path GridKit::emt_path GridKit::testing)

add_test(NAME EMTPathTest COMMAND $<TARGET_FILE:test_emt_path>)

add_executable(test_emt_tower runTowerTests.cpp)
target_link_libraries(
  test_emt_tower GridKit::emt_tower GridKit::testing)

add_test(NAME EMTTowerTest COMMAND $<TARGET_FILE:test_emt_tower>)

add_executable(test_emt_frequency_response_parser runFrequencyResponseParserTests.cpp)
target_link_libraries(
  test_emt_frequency_response_parser GridKit::emt_overhead GridKit::testing)
target_include_directories(
  test_emt_frequency_response_parser
  PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/nlohmann-json/include)

add_test(NAME EMTFrequencyResponseParserTest
         COMMAND $<TARGET_FILE:test_emt_frequency_response_parser>)

add_executable(test_emt_dense_expression runDenseExpressionTests.cpp)
target_link_libraries(
  test_emt_dense_expression GridKit::emt_core GridKit::testing)

add_test(NAME EMTDenseExpressionTest COMMAND $<TARGET_FILE:test_emt_dense_expression>)

add_executable(test_emt_element_jacobian runElementJacobianTests.cpp)
target_link_libraries(
  test_emt_element_jacobian GridKit::emt_core GridKit::testing)

add_test(NAME EMTElementJacobianTest COMMAND $<TARGET_FILE:test_emt_element_jacobian>)

add_executable(test_emt_skin_effect runSkinEffectTests.cpp)
target_link_libraries(
  test_emt_skin_effect GridKit::emt_skin_effect GridKit::testing)

add_test(NAME EMTSkinEffectTest COMMAND $<TARGET_FILE:test_emt_skin_effect>)
