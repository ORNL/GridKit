#

add_executable(test_ida runIdaTests.cpp)
target_link_libraries(
  test_ida GridKit::solvers_dyn GridKit::testing)

# add_test(NAME IDATest COMMAND $<TARGET_FILE:test_ida>)

install(TARGETS test_ida RUNTIME DESTINATION bin)

if(GRIDKIT_ENABLE_SUNDIALS_SPARSE AND TARGET SUNDIALS::arkode)
  add_executable(test_arkode runArkodeTests.cpp)
  target_link_libraries(test_arkode GridKit::solvers_dyn GridKit::testing)

  add_test(NAME ArkodeTest COMMAND $<TARGET_FILE:test_arkode>)

  install(TARGETS test_arkode RUNTIME DESTINATION bin)
endif()
