#

if(TARGET SUNDIALS::idas)
  add_executable(test_ida runIdaTests.cpp)
  target_link_libraries(test_ida GridKit::solvers_dyn GridKit::testing)

  add_test(NAME IDATest COMMAND $<TARGET_FILE:test_ida>)

  install(TARGETS test_ida RUNTIME DESTINATION bin)
endif()

if(GRIDKIT_ENABLE_RESOLVE)
  add_executable(test_implicit_trapezoidal runImplicitTrapezoidalTests.cpp)
  target_link_libraries(
    test_implicit_trapezoidal
    GridKit::implicit_trapezoidal
    GridKit::rms_norm
    GridKit::testing
    GridKit::resolve_system_solver
    ReSolve::ReSolve)
  add_test(NAME ImplicitTrapezoidalTest COMMAND $<TARGET_FILE:test_implicit_trapezoidal>)

  install(TARGETS test_implicit_trapezoidal RUNTIME DESTINATION bin)

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
