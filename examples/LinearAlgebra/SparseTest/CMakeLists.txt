

add_executable(spmattest SparseTest.cpp)
target_link_libraries(spmattest GridKit::SparseMatrix)
                                
add_test(NAME SparseMatrixTest    COMMAND $<TARGET_FILE:spmattest>)
install(TARGETS spmattest RUNTIME DESTINATION bin)
