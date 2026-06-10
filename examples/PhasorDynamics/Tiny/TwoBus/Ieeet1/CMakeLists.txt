gridkit_example_current_install_path(_install_path)

add_executable(TwoBusIeeet1 TwoBusIeeet1.cpp)
target_link_libraries(
  TwoBusIeeet1 GridKit::phasor_dynamics_systemmodel GridKit::solvers_dyn)
install(TARGETS TwoBusIeeet1 RUNTIME DESTINATION ${_install_path})

add_executable(TwoBusIeeet1Json TwoBusIeeet1Json.cpp)
target_link_libraries(
  TwoBusIeeet1Json GridKit::phasor_dynamics_systemmodel GridKit::solvers_dyn)
install(TARGETS TwoBusIeeet1Json RUNTIME DESTINATION ${_install_path})
gridkit_example_add_file(TwoBusIeeet1.json)

add_test(NAME TwoBusIeeet1 COMMAND TwoBusIeeet1)
add_test(
  NAME TwoBusIeeet1Json COMMAND TwoBusIeeet1Json ${CMAKE_CURRENT_BINARY_DIR}/TwoBusIeeet1.json)
add_test(
  NAME TwoBusIeeet1Json_no_arg
  COMMAND TwoBusIeeet1Json
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
