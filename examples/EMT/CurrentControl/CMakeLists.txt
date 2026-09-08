foreach(
  _file
  GFL.solver.json
  GFM.solver.json
  README.md
  run.py
  validate.py
  plot.py)
  configure_file(${_file} ${_file} COPYONLY)
  install(FILES ${_file} DESTINATION share/gridkit/examples/EMT/CurrentControl)
endforeach()

foreach(
  _file
  GFL.case.json
  GFM.case.json
  GFL.state.json
  GFM.state.json
  README.md)
  configure_file(${CMAKE_SOURCE_DIR}/cases/EMT/CurrentControl/${_file}
                 ${CMAKE_BINARY_DIR}/cases/EMT/CurrentControl/${_file} COPYONLY)
  install(FILES ${CMAKE_SOURCE_DIR}/cases/EMT/CurrentControl/${_file}
          DESTINATION share/gridkit/cases/EMT/CurrentControl)
endforeach()

find_package(
  Python3
  COMPONENTS Interpreter
  QUIET)
if(Python3_Interpreter_FOUND
   AND GRIDKIT_ENABLE_ENZYME
   AND GRIDKIT_ENABLE_SUNDIALS_SPARSE)
  add_test(
    NAME EMTGflPowerControl
    COMMAND
      ${Python3_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/validate.py --exe
      $<TARGET_FILE:EMTDynamicSimulation>)
  set_tests_properties(EMTGflPowerControl PROPERTIES TIMEOUT 300)
endif()
