add_executable(EnzymeLibVectorCheck EnzymeVector.cpp VectorModel.cpp)
target_link_libraries(EnzymeLibVectorCheck
    ClangEnzymeFlags
    GridKit::DenseMatrix
    GridKit::Utilities)

add_test(NAME "EnzymeLibVectorCheck" COMMAND ${CMAKE_CURRENT_BINARY_DIR}/EnzymeLibVectorCheck)
