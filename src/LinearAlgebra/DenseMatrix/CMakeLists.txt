add_library(DenseMatrix INTERFACE)
include_directories(DenseMatrix INTERFACE ${CMAKE_CURRENT_LIST_DIR})

add_library(GRIDKIT::DenseMatrix ALIAS DenseMatrix)
