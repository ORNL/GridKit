#

if(TARGET SUNDIALS::idas)
  add_executable(test_ida runIdaTests.cpp)
  target_link_libraries(test_ida GridKit::solvers_dyn GridKit::testing)

  add_test(NAME IDATest COMMAND $<TARGET_FILE:test_ida>)

  add_executable(test_partitioned runPartitionedTests.cpp)
  target_link_libraries(test_partitioned GridKit::solvers_dyn GridKit::testing)
  add_test(NAME PartitionedTest COMMAND $<TARGET_FILE:test_partitioned>)

  install(TARGETS test_ida test_partitioned RUNTIME DESTINATION bin)
endif()

if(GRIDKIT_ENABLE_RESOLVE)
  add_executable(test_rosenbrock runRosenbrockTests.cpp)
  target_link_libraries(
    test_rosenbrock
    GridKit::rosenbrock
    GridKit::testing
    GridKit::fixed_step
    GridKit::resolve_system_solver
    ReSolve::ReSolve)
  add_test(NAME RosenbrockTest COMMAND $<TARGET_FILE:test_rosenbrock>)

  install(TARGETS test_rosenbrock RUNTIME DESTINATION bin)
endif()
