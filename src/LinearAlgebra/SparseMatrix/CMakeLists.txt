
add_library(SparseMatrix INTERFACE)
include_directories(SparseMatrix INTERFACE ${CMAKE_CURRENT_LIST_DIR})

add_library(GRIDKIT::SparseMatrix ALIAS SparseMatrix)

install(
  FILES
    COO_Matrix.hpp
  DESTINATION
    include/GridKit/LinearAlgebra/SparseMatrix)
