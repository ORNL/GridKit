foreach(
  _file
  Circuit.case.json
  Hybrid.case.json
  Hybrid.solver.json
  README.md
  validate.py)
  configure_file(${CMAKE_SOURCE_DIR}/cases/EMT/FrequencyDependent/${_file}
                 ${CMAKE_BINARY_DIR}/cases/EMT/FrequencyDependent/${_file} COPYONLY)
  install(FILES ${CMAKE_SOURCE_DIR}/cases/EMT/FrequencyDependent/${_file}
          DESTINATION share/gridkit/cases/EMT/FrequencyDependent)
endforeach()

find_package(
  Python3
  COMPONENTS Interpreter
  QUIET)
if(Python3_Interpreter_FOUND
   AND GRIDKIT_ENABLE_ENZYME
   AND GRIDKIT_ENABLE_SUNDIALS_SPARSE)
  add_test(
    NAME EMTFrequencyDependent
    COMMAND
      ${Python3_EXECUTABLE} ${CMAKE_BINARY_DIR}/cases/EMT/FrequencyDependent/validate.py --exe
      $<TARGET_FILE:EMTDynamicSimulation>)
  set_tests_properties(EMTFrequencyDependent PROPERTIES TIMEOUT 300)
endif()
