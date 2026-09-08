foreach(
  _file
  Hawaii.case.json
  Hawaii.state.json
  Hawaii.solver.json
  conversion.json
  README.md
  convert.py
  validate.py)
  configure_file(${CMAKE_SOURCE_DIR}/cases/EMT/Hawaii/${_file}
                 ${CMAKE_BINARY_DIR}/cases/EMT/Hawaii/${_file} COPYONLY)
  install(FILES ${CMAKE_SOURCE_DIR}/cases/EMT/Hawaii/${_file}
          DESTINATION share/gridkit/cases/EMT/Hawaii)
endforeach()

find_package(
  Python3
  COMPONENTS Interpreter
  QUIET)
if(Python3_Interpreter_FOUND
   AND GRIDKIT_ENABLE_ENZYME
   AND GRIDKIT_ENABLE_SUNDIALS_SPARSE)
  add_test(
    NAME EMTHawaiiCase
    COMMAND
      ${Python3_EXECUTABLE} ${CMAKE_BINARY_DIR}/cases/EMT/Hawaii/validate.py --exe
      $<TARGET_FILE:EMTDynamicSimulation>)
  set_tests_properties(EMTHawaiiCase PROPERTIES TIMEOUT 1800)
endif()
