foreach(
  _file
  DCLink.case.json
  DCLink.state.json
  DCLink.solver.json
  README.md
  validate.py)
  configure_file(${CMAKE_SOURCE_DIR}/cases/EMT/DCLink/${_file}
                 ${CMAKE_BINARY_DIR}/cases/EMT/DCLink/${_file} COPYONLY)
  install(FILES ${CMAKE_SOURCE_DIR}/cases/EMT/DCLink/${_file}
          DESTINATION share/gridkit/cases/EMT/DCLink)
endforeach()

find_package(
  Python3
  COMPONENTS Interpreter
  QUIET)
if(Python3_Interpreter_FOUND
   AND GRIDKIT_ENABLE_ENZYME
   AND GRIDKIT_ENABLE_SUNDIALS_SPARSE)
  add_test(
    NAME EMTDcLinkTransient
    COMMAND
      ${Python3_EXECUTABLE} ${CMAKE_BINARY_DIR}/cases/EMT/DCLink/validate.py --exe
      $<TARGET_FILE:EMTDynamicSimulation>)
  set_tests_properties(EMTDcLinkTransient PROPERTIES TIMEOUT 120)
endif()
