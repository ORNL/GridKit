

add_executable(spmattest SparseTest.cpp)
target_link_libraries(spmattest GRIDKIT::SparseMatrix)
                                
add_test(NAME SparseMatrixTest    COMMAND $<TARGET_FILE:spmattest>)
install(TARGETS spmattest RUNTIME DESTINATION bin)
