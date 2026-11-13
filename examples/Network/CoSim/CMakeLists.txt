gridkit_example_add_file(TwoBusCoSimClient.case.json)
gridkit_example_add_file(TwoBusCoSimClient.solo.case.json)
gridkit_example_add_file(TwoBusCoSimClient.solver.json)
gridkit_example_add_file(TwoBusCoSimServer.case.json)
gridkit_example_add_file(TwoBusCoSimServer.solo.case.json)
gridkit_example_add_file(TwoBusCoSimServer.solver.json)
gridkit_example_add_file(ThreeBusCoSimClient.case.json)
gridkit_example_add_file(ThreeBusCoSimServer.case.json)
gridkit_example_add_file(IEEE39Generator.case.json)
gridkit_example_add_file(IEEE39Grid.case.json)
gridkit_example_add_file(run_cosim.sh)

add_executable(CoSimClient CoSimClient.cpp)
target_link_libraries(
  CoSimClient
  PUBLIC GridKit::utilities_cli_args
         GridKit::phasor_dynamics_systemmodel
         GridKit::solvers_dyn
         cppzmq)

add_executable(CoSimServer CoSimServer.cpp)
target_link_libraries(
  CoSimServer
  PUBLIC GridKit::utilities_cli_args
         GridKit::phasor_dynamics_systemmodel
         GridKit::solvers_dyn
         cppzmq)

foreach(_target CoSimClient CoSimServer)
  target_include_directories(${_target}
                             PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/nlohmann-json/single_include)
endforeach()

add_test(
  NAME CoSimExample
  COMMAND bash run_cosim.sh
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

set_tests_properties(CoSimExample PROPERTIES TIMEOUT 30)

add_test(
  NAME IEEE39CoSimExample
  COMMAND bash run_cosim.sh ieee39
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(
  IEEE39CoSimExample
  PROPERTIES TIMEOUT
             30
             RESOURCE_LOCK
             cosim_port_5556)
set_tests_properties(CoSimExample PROPERTIES RESOURCE_LOCK cosim_port_5556)

add_test(
  NAME CoSimSignalsExample
  COMMAND bash run_cosim.sh three-bus signals
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(
  CoSimSignalsExample
  PROPERTIES TIMEOUT
             30
             RESOURCE_LOCK
             cosim_port_5556)
