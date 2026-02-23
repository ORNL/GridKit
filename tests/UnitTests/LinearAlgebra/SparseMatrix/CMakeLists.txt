add_executable(test_sparse_csr runSparseCsrTests.cpp)
target_link_libraries(test_sparse_csr GridKit::sparse_matrix GridKit::testing)

add_test(NAME SparseTest COMMAND $<TARGET_FILE:test_sparse_csr>)

add_executable(test_sparse_coo runSparseCooTests.cpp)
target_link_libraries(test_sparse_coo GridKit::sparse_matrix GridKit::testing)

add_test(NAME SparseCooTest COMMAND $<TARGET_FILE:test_sparse_coo>)

install(TARGETS test_sparse_csr test_sparse_coo
        RUNTIME DESTINATION bin)
