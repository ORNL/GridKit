foreach(
  _file
  RL.case.json
  Steady.state.json
  Perturbed.state.json
  Default.state.json
  Steady.solver.json
  Perturbed.solver.json
  Default.solver.json
  README.md
  validate.py)
  configure_file(${CMAKE_SOURCE_DIR}/cases/EMT/CurrentInitialization/${_file}
                 ${CMAKE_BINARY_DIR}/cases/EMT/CurrentInitialization/${_file} COPYONLY)
  install(FILES ${CMAKE_SOURCE_DIR}/cases/EMT/CurrentInitialization/${_file}
          DESTINATION share/gridkit/cases/EMT/CurrentInitialization)
endforeach()

find_package(
  Python3
  COMPONENTS Interpreter
  QUIET)
if(Python3_Interpreter_FOUND
   AND GRIDKIT_ENABLE_ENZYME
   AND GRIDKIT_ENABLE_SUNDIALS_SPARSE)
  add_test(
    NAME EMTCurrentInitialization
    COMMAND
      ${Python3_EXECUTABLE} ${CMAKE_BINARY_DIR}/cases/EMT/CurrentInitialization/validate.py --exe
      $<TARGET_FILE:EMTDynamicSimulation>)
  set_tests_properties(EMTCurrentInitialization PROPERTIES TIMEOUT 120)
endif()
