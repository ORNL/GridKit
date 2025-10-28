
add_library(SparseMatrix CsrMatrix.cpp)
include_directories(SparseMatrix ${CMAKE_CURRENT_LIST_DIR})

add_library(GRIDKIT::SparseMatrix ALIAS SparseMatrix)
