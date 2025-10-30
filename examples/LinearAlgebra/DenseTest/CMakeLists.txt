
add_executable(densemattest DenseTest.cpp)
target_link_libraries(densemattest GridKit::DenseMatrix)
                                
add_test(NAME DenseMatrixTest COMMAND $<TARGET_FILE:densemattest>)
install(TARGETS densemattest RUNTIME DESTINATION bin)
