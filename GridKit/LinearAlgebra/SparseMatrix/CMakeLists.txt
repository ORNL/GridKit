add_library(CsrMatrix INTERFACE)
target_include_directories(CsrMatrix INTERFACE
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>)

add_library(GridKit::CsrMatrix ALIAS CsrMatrix)

install(
  FILES
    COO_Matrix.hpp
    CsrMatrix.hpp
  DESTINATION
    include/GridKit/LinearAlgebra/SparseMatrix)
