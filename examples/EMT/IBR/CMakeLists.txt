foreach(
  _file
  01_Baseline.solver.json
  02_LoadStep.solver.json
  03_LoadPulse.solver.json
  04_FaultClearing.solver.json
  05_TieReclosing.solver.json
  06_FaultClearingSwitching.solver.json
  run.py
  plot.py
  README.md)
  gridkit_example_add_file(${_file})
endforeach()

# Preserve the solver files' relative case paths in build and install trees.
foreach(
  _file
  TenBus.case.json
  TenBus.state.json
  README.md
  build_case.py)
  configure_file(${CMAKE_SOURCE_DIR}/cases/EMT/IBR/${_file}
                 ${CMAKE_BINARY_DIR}/cases/EMT/IBR/${_file} COPYONLY)
  install(FILES ${CMAKE_SOURCE_DIR}/cases/EMT/IBR/${_file} DESTINATION share/gridkit/cases/EMT/IBR)
endforeach()

find_package(
  Python3
  COMPONENTS Interpreter
  QUIET)
if(Python3_Interpreter_FOUND
   AND GRIDKIT_ENABLE_ENZYME
   AND GRIDKIT_ENABLE_SUNDIALS_SPARSE)
  add_test(
    NAME EMTTenBusIBR
    COMMAND
      ${Python3_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/run.py --exe
      $<TARGET_FILE:EMTDynamicSimulation> --results ${CMAKE_CURRENT_BINARY_DIR}/smoke --smoke)
  set_tests_properties(EMTTenBusIBR PROPERTIES TIMEOUT 180)
endif()
