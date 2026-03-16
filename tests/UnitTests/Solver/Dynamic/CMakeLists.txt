#

add_executable(test_ida runIdaTests.cpp)
target_link_libraries(
  test_ida GridKit::solvers_dyn GridKit::testing)

add_test(NAME IDATest COMMAND $<TARGET_FILE:test_ida>)

add_executable(test_ros runRosenbrockTests.cpp)
target_link_libraries(test_ros GridKit::rosenbrock GridKit::testing)
add_test(NAME ROSTest COMMAND $<TARGET_FILE:test_ros>)
                                

install(TARGETS test_ida RUNTIME DESTINATION bin)
