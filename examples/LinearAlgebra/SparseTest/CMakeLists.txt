
add_executable(spmattest SparseTest.cpp)
target_link_libraries(
  spmattest GridKit::sparse_matrix)

add_test(NAME SparseMatrixTest COMMAND spmattest)
install(TARGETS spmattest RUNTIME DESTINATION bin)
