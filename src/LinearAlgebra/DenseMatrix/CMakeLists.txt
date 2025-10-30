
add_library(DenseMatrix INTERFACE)
include_directories(DenseMatrix INTERFACE ${CMAKE_CURRENT_LIST_DIR})

add_library(GridKit::DenseMatrix ALIAS DenseMatrix)

install(
  FILES
    DenseMatrix.hpp
  DESTINATION
    include/GridKit/LinearAlgebra/DenseMatrix)
