add_library(DenseMatrix INTERFACE)
target_include_directories(
  DenseMatrix
  INTERFACE ${CMAKE_SOURCE_DIR})

add_library(GridKit::DenseMatrix ALIAS DenseMatrix)

install(FILES DenseMatrix.hpp DESTINATION include/GridKit/LinearAlgebra/DenseMatrix)
