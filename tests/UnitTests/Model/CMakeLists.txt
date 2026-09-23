add_executable(test_state_data runStateDataTests.cpp)
target_link_libraries(
  test_state_data
  PRIVATE GridKit::model_state_data GridKit::testing)

add_test(NAME ModelStateDataTest COMMAND test_state_data)

install(TARGETS test_state_data)
