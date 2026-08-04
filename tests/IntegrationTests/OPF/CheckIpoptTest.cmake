if(NOT DEFINED TEST_EXECUTABLE)
  message(FATAL_ERROR "TEST_EXECUTABLE is required")
endif()

execute_process(
  COMMAND "${TEST_EXECUTABLE}"
  RESULT_VARIABLE test_result
  OUTPUT_VARIABLE test_output
  ERROR_VARIABLE test_error)

message("${test_output}")
if(test_error)
  message("${test_error}")
endif()

if(NOT
   test_result
   EQUAL
   0)
  message(FATAL_ERROR "OPF Ipopt test exited with status ${test_result}")
endif()

if(NOT
   test_output
   MATCHES
   "Starting derivative checker for second derivatives")
  message(FATAL_ERROR "Ipopt second-order derivative checker did not run")
endif()

if(NOT
   test_output
   MATCHES
   "No errors detected by derivative checker\\.")
  message(FATAL_ERROR "Ipopt derivative checker did not pass")
endif()
