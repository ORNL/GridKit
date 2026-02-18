add_executable(test_sparse runSparseTests.cpp)
target_link_libraries(test_sparse GridKit::sparse)

add_test(NAME SparseTest COMMAND $<TARGET_FILE:test_sparse>)

add_executable(test_sparse_coo runSparseCooTests.cpp)
target_link_libraries(test_sparse_coo GridKit::sparse)

add_test(NAME SparseCooTest COMMAND $<TARGET_FILE:test_sparse_coo>)

install(TARGETS test_sparse test_sparse_coo
        RUNTIME DESTINATION bin)
