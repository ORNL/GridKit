add_executable(EnzymeLibScalarCheck EnzymeScalar.cpp ScalarModel.cpp)
target_link_libraries(EnzymeLibScalarCheck ClangEnzymeFlags)

add_test(NAME "EnzymeLibScalarCheck" COMMAND ${CMAKE_CURRENT_BINARY_DIR}/EnzymeLibScalarCheck)
