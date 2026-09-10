# [[
#  Author(s):
#    - Slaven Peles <peless@ornl.gov>
#]]

add_subdirectory(Dynamic)

if(TARGET GridKit::solvers_steady)
  add_executable(test_kinsol runKinsolTests.cpp)
  target_link_libraries(
    test_kinsol
    GridKit::solvers_steady
    GridKit::testing)
  add_test(NAME KinsolTest COMMAND test_kinsol)
endif()
