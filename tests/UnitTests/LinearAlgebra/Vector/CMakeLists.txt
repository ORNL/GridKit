add_executable(test_vector runVectorTests.cpp)
target_link_libraries(
  test_vector GridKit::dense_vector GridKit::testing)

add_test(NAME VectorTest COMMAND $<TARGET_FILE:test_vector>)

install(TARGETS test_vector RUNTIME DESTINATION bin)
