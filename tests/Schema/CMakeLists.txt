# Line-input schema validation; runs whenever a Python interpreter with
# jsonschema is available so the schema and the example corpus cannot
# drift silently.

find_package(Python3 COMPONENTS Interpreter)
if(Python3_Interpreter_FOUND)
  execute_process(
    COMMAND ${Python3_EXECUTABLE} -c "import jsonschema"
    RESULT_VARIABLE _jsonschema_missing
    OUTPUT_QUIET ERROR_QUIET)
  if(_jsonschema_missing EQUAL 0)
    add_test(
      NAME LineSchemaValidation
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/validate_line_inputs.py)
  else()
    message(STATUS "jsonschema not found; skipping LineSchemaValidation")
  endif()
endif()
