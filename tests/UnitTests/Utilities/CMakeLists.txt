add_executable(test_case_format runCaseFormatTests.cpp)
target_link_libraries(test_case_format PRIVATE
    GRIDKIT::phasor_dynamics_components)
target_include_directories(test_case_format
    PRIVATE
    ${GRIDKIT_THIRD_PARTY_DIR}/nlohmann-json/single_include
    ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include)

add_executable(test_logger runLoggerTests.cpp)
target_link_libraries(test_logger PRIVATE GRIDKIT::utilities_logger)
add_executable(test_cli_options runCliOptionsTests.cpp)
target_link_libraries(test_cli_options PRIVATE GRIDKIT::utilities_cli_options)

add_test(NAME UtilitiesCaseFormatTest COMMAND test_case_format)
add_test(NAME UtilitiesLoggerTest COMMAND test_logger)
add_test(NAME UtilitiesCliOptionsTest COMMAND test_cli_options)

install(TARGETS test_case_format)
install(TARGETS test_logger)
install(TARGETS test_cli_options)
