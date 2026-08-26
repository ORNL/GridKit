gridkit_example_add_file(TwoBusCoSimClient.case.json)
gridkit_example_add_file(TwoBusCoSimServer.case.json)
gridkit_example_add_file(ThreeBusCoSimClient.case.json)
gridkit_example_add_file(ThreeBusCoSimServer.case.json)
gridkit_example_add_file(run_cosim.sh)

add_executable(CoSimClient CoSimClient.cpp)
target_link_libraries(
  CoSimClient
  PUBLIC GridKit::utilities_cli_args GridKit::phasor_dynamics_systemmodel GridKit::solvers_dyn ZMQ)

add_executable(CoSimServer CoSimServer.cpp)
target_link_libraries(
  CoSimServer
  PUBLIC GridKit::utilities_cli_args GridKit::phasor_dynamics_systemmodel GridKit::solvers_dyn ZMQ)

add_test(
  NAME CoSimExample
  COMMAND bash run_cosim.sh
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
