# 

add_executable(test_csr runCsrTests.cpp)
target_include_directories(test_csr PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>)
                                

add_test(NAME IDATest COMMAND $<TARGET_FILE:test_csr>)

install(TARGETS test_csr 
        RUNTIME DESTINATION bin)
