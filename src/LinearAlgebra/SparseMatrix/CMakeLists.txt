
add_library(CsrMatrix CsrMatrix.cpp)

install(
  FILES
    COO_Matrix.hpp
  DESTINATION
    include/GridKit/LinearAlgebra/SparseMatrix)
	
add_library(GridKit::CsrMatrix ALIAS CsrMatrix)
