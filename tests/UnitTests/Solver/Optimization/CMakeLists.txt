if(TARGET solvers_opt)
  add_executable(test_ipopt_evaluator runIpoptEvaluatorTests.cpp)
  target_link_libraries(test_ipopt_evaluator GridKit::solvers_opt)

  add_test(NAME IpoptEvaluatorTest COMMAND $<TARGET_FILE:test_ipopt_evaluator>)
endif()
