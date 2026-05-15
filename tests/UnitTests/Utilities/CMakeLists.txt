add_executable(test_case_format runCaseFormatTests.cpp)
target_link_libraries(test_case_format PRIVATE
    GridKit::phasor_dynamics_components
    GridKit::testing)
target_include_directories(test_case_format
    PRIVATE
    ${GRIDKIT_THIRD_PARTY_DIR}/nlohmann-json/single_include
    ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include)

add_executable(test_logger runLoggerTests.cpp)
target_link_libraries(test_logger PRIVATE
    GridKit::utilities_logger
    GridKit::testing)
add_executable(test_cli_args runCliArgsTests.cpp)
target_link_libraries(test_cli_args PRIVATE
    GridKit::utilities_cli_args
    GridKit::utilities_logger
    GridKit::testing)

add_test(NAME UtilitiesCaseFormatTest COMMAND test_case_format)
add_test(NAME UtilitiesLoggerTest COMMAND test_logger)
add_test(NAME UtilitiesCliArgsTest COMMAND test_cli_args)

install(TARGETS test_case_format)
install(TARGETS test_logger)
install(TARGETS test_cli_args)
