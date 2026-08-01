add_executable(test_opf_ipopt runIpoptSystemTests.cpp)
target_link_libraries(
  test_opf_ipopt
  PRIVATE GridKit::model_opf GridKit::solvers_opt)

add_test(NAME OPFIpoptSystemTest COMMAND test_opf_ipopt)

install(TARGETS test_opf_ipopt)
