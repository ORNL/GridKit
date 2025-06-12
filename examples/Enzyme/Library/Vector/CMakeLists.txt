add_executable(EnzymeLibVectorCheck EnzymeVector.cpp VectorModel.cpp)
target_link_libraries(EnzymeLibVectorCheck ClangEnzymeFlags GRIDKIT::DenseMatrix GRIDKIT::Utilities)

add_test(NAME "EnzymeLibVectorCheck" COMMAND ${CMAKE_CURRENT_BINARY_DIR}/EnzymeLibVectorCheck)
