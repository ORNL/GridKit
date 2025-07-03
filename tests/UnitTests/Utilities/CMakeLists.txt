add_executable(test_case_format runCaseFormatTests.cpp)

add_test(NAME UtilitiesCaseFormatTest COMMAND $<TARGET_FILE:test_case_format>)

install(TARGETS test_case_format)
