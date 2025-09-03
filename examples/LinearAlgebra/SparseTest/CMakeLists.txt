

add_executable(spmattest SparseTest.cpp)
                                
add_test(NAME SparseMatrixTest COMMAND spmattest)
install(TARGETS spmattest RUNTIME DESTINATION bin)
