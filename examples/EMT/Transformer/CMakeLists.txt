foreach(
  _file
  Transformer.case.json
  Steady.state.json
  Energization.state.json
  Steady.solver.json
  Energization.solver.json
  README.md
  steady_state.py
  validate.py)
  configure_file(${CMAKE_SOURCE_DIR}/cases/EMT/Transformer/${_file}
                 ${CMAKE_BINARY_DIR}/cases/EMT/Transformer/${_file} COPYONLY)
  install(FILES ${CMAKE_SOURCE_DIR}/cases/EMT/Transformer/${_file}
          DESTINATION share/gridkit/cases/EMT/Transformer)
endforeach()

find_package(
  Python3
  COMPONENTS Interpreter
  QUIET)
if(Python3_Interpreter_FOUND
   AND GRIDKIT_ENABLE_ENZYME
   AND GRIDKIT_ENABLE_SUNDIALS_SPARSE)
  add_test(
    NAME EMTTransformerCase
    COMMAND
      ${Python3_EXECUTABLE} ${CMAKE_BINARY_DIR}/cases/EMT/Transformer/validate.py --exe
      $<TARGET_FILE:EMTDynamicSimulation>)
  set_tests_properties(EMTTransformerCase PROPERTIES TIMEOUT 300)
endif()
