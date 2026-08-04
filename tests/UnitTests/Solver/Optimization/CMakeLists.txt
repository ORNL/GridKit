add_executable(test_optimization_ipopt runIpoptTests.cpp)
target_link_libraries(
  test_optimization_ipopt
  PRIVATE test_optimization_model GridKit::solvers_opt_ipopt GridKit::testing)
set_property(TARGET test_optimization_ipopt PROPERTY LINKER_LANGUAGE Fortran)

add_test(
  NAME OptimizationIpoptExactTest
  COMMAND
    ${CMAKE_COMMAND} -DTEST_EXECUTABLE=$<TARGET_FILE:test_optimization_ipopt> -P
    ${CMAKE_CURRENT_SOURCE_DIR}/CheckIpoptTest.cmake)

install(TARGETS test_optimization_ipopt RUNTIME DESTINATION bin)
