if(GRIDKIT_ENABLE_ZMQ)
  add_executable(test_zmq runZMQTests.cpp)
  target_link_libraries(test_zmq PRIVATE GridKit::testing cppzmq)

  add_test(NAME PhasorDynamicsZMQTest COMMAND test_zmq)
  install(TARGETS test_zmq RUNTIME DESTINATION bin)
endif()
