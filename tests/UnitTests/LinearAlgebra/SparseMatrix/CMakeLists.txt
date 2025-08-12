# 

add_executable(test_csr runCSRTests.cpp)
target_link_libraries(test_csr GRIDKIT::SparseMatrix)
                                

add_test(NAME IDATest COMMAND $<TARGET_FILE:test_csr>)

install(TARGETS test_csr 
        RUNTIME DESTINATION bin)
