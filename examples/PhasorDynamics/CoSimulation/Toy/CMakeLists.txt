foreach(_case TwoBusCoSimClient TwoBusCoSimServer
              ThreeBusCoSimClient ThreeBusCoSimServer)
  gridkit_example_add_file(${_case}.case.json)
  gridkit_example_add_file(${_case}.solver.json)
endforeach()
