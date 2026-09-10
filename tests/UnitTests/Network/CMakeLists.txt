if(GRIDKIT_ENABLE_ZMQ)
  add_executable(test_zmq runZMQTests.cpp)
  target_link_libraries(test_zmq PRIVATE GridKit::testing cppzmq)

  add_test(NAME PhasorDynamicsZMQTest COMMAND test_zmq)
  install(TARGETS test_zmq RUNTIME DESTINATION bin)

  add_executable(test_cosim runCoSimTests.cpp)
  target_link_libraries(
    test_cosim
    GridKit::co_sim
    GridKit::testing
    cppzmq)
  target_include_directories(
    test_cosim
    PRIVATE
    ${GRIDKIT_THIRD_PARTY_DIR}/nlohmann-json/include)

  add_test(NAME CoSimulationChannelTest COMMAND test_cosim)
  set_tests_properties(CoSimulationChannelTest PROPERTIES TIMEOUT 30)
endif()
