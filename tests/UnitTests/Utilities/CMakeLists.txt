add_executable(test_case_format runCaseFormatTests.cpp)
target_include_directories(test_case_format
    PRIVATE
    ${GRIDKIT_THIRD_PARTY_DIR}/nlohmann-json/single_include
    ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include)

add_executable(test_logger runLoggerTests.cpp)
target_link_libraries(test_logger PRIVATE GRIDKIT::utilities_logger)

add_test(NAME UtilitiesCaseFormatTest COMMAND test_case_format)
add_test(NAME UtilitiesLoggerTest COMMAND test_logger)

install(TARGETS test_case_format)
install(TARGETS test_logger)
