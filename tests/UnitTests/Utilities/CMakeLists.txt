add_executable(test_case_format runCaseFormatTests.cpp)
target_link_libraries(test_case_format PRIVATE nlohmann_json::nlohmann_json magic_enum::magic_enum)

add_test(NAME UtilitiesCaseFormatTest COMMAND $<TARGET_FILE:test_case_format>)

install(TARGETS test_case_format)
