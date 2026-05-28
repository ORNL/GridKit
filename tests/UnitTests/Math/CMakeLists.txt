add_executable(test_math_smoothness_indicator runSmoothnessIndicatorTests.cpp)
target_link_libraries(
  test_math_smoothness_indicator
  PRIVATE GridKit::definitions GridKit::testing)

add_test(NAME MathSmoothnessIndicatorTest COMMAND test_math_smoothness_indicator)

install(TARGETS test_math_smoothness_indicator)
