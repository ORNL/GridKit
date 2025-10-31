add_executable(test_sparse runSparseTests.cpp)
target_link_libraries(test_sparse GRIDKIT::CsrMatrix)
                                
add_test(NAME SparseTest COMMAND $<TARGET_FILE:test_sparse>)


install(TARGETS test_sparse 
        RUNTIME DESTINATION bin)
