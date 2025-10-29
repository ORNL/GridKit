
add_library(SparseMatrix CsrMatrix.cpp)

add_library(GridKit::SparseMatrix ALIAS SparseMatrix)

install(
  FILES
    COO_Matrix.hpp
  DESTINATION
    include/GridKit/LinearAlgebra/SparseMatrix)
