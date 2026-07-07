#

if(TARGET SUNDIALS::idas)
  add_executable(test_ida runIdaTests.cpp)
  target_link_libraries(test_ida GridKit::solvers_dyn GridKit::testing)

  add_test(NAME IDATest COMMAND $<TARGET_FILE:test_ida>)

  install(TARGETS test_ida RUNTIME DESTINATION bin)
endif()

if(GRIDKIT_ENABLE_RESOLVE)
  add_executable(test_ros runRosenbrockTests.cpp)
  target_link_libraries(test_ros GridKit::rosenbrock GridKit::testing)
  add_test(NAME ROSTest COMMAND $<TARGET_FILE:test_ros>)

  install(TARGETS test_ros RUNTIME DESTINATION bin)
endif()
