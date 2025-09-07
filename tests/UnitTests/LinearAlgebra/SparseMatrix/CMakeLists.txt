# 

add_executable(test_csr runCsrTests.cpp)
target_include_directories(test_csr PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>)
                                

add_test(NAME IDATest COMMAND $<TARGET_FILE:test_csr>)

add_executable(test_jac runCsrJacobianTests.cpp)
target_link_libraries(test_jac
        GridKit::power_elec_resistor)
add_test(NAME CSRJacTest COMMAND $<TARGET_FILE:test_jac>)

install(TARGETS test_csr 
        RUNTIME DESTINATION bin)
