enzyme_add_executable(
  NAME EnzymeLibVectorCheck 
  SOURCES EnzymeVector.cpp VectorModel.cpp
  LINK_LIBRARIES GRIDKIT::DenseMatrix
)

add_test(NAME "EnzymeLibVectorCheck" COMMAND ${CMAKE_CURRENT_BINARY_DIR}/EnzymeLibVectorCheck)
