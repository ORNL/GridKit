
add_library(SparseMatrix CsrMatrix.cpp)
include_directories(SparseMatrix ${CMAKE_CURRENT_LIST_DIR})

add_library(GridKit::SparseMatrix ALIAS SparseMatrix)

install(
  FILES
    COO_Matrix.hpp
  DESTINATION
    include/GridKit/LinearAlgebra/SparseMatrix)
