
add_library(SparseMatrix INTERFACE)
target_include_directories(SparseMatrix INTERFACE ${CMAKE_SOURCE_DIR})

add_library(GridKit::SparseMatrix ALIAS SparseMatrix)

install(
  FILES
    COO_Matrix.hpp
    CsrMatrix.hpp
    CsrPermutedSum.hpp
  DESTINATION
    include/GridKit/LinearAlgebra/SparseMatrix)