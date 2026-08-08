foreach(validation_case IN ITEMS ACTIVSg200 ACTIVSg500 WECC240)
  set(validation_work_dir ${CMAKE_CURRENT_BINARY_DIR}/${validation_case})
  file(MAKE_DIRECTORY ${validation_work_dir})

  add_test(
    NAME ${validation_case}_validation
    COMMAND
      DynamicSimulation
      ${CMAKE_CURRENT_SOURCE_DIR}/${validation_case}/${validation_case}.solver.json
    WORKING_DIRECTORY ${validation_work_dir})
endforeach()
