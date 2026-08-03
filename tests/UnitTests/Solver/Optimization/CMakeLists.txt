# Optimization family unit tests; the fitter needs the Ipopt-backed
# quadratic-program stack.

if(TARGET GridKit::vector_fitting)
  add_executable(test_vector_fitting runVectorFittingTests.cpp)
  target_link_libraries(test_vector_fitting GridKit::vector_fitting GridKit::testing)

  add_test(NAME VectorFittingTest COMMAND $<TARGET_FILE:test_vector_fitting>)

  install(TARGETS test_vector_fitting RUNTIME DESTINATION bin)
endif()
