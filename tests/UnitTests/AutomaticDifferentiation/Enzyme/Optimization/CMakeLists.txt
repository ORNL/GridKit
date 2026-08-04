add_library(test_optimization_model STATIC ModelEnzyme.cpp)
target_include_directories(
  test_optimization_model
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(
  test_optimization_model
  PUBLIC GridKit::optimization_core
  PRIVATE GridKit::optimization_enzyme)
target_compile_features(test_optimization_model PUBLIC cxx_std_20)

add_executable(test_optimization_enzyme runOptimizationTests.cpp)
target_link_libraries(
  test_optimization_enzyme
  PRIVATE test_optimization_model GridKit::testing)

add_test(NAME OptimizationEnzymeDerivativeTest COMMAND $<TARGET_FILE:test_optimization_enzyme>)

install(TARGETS test_optimization_enzyme RUNTIME DESTINATION bin)
