
add_library(SparseMatrix INTERFACE)
target_include_directories(SparseMatrix INTERFACE ${CMAKE_SOURCE_DIR})

add_library(GRIDKIT::SparseMatrix ALIAS SparseMatrix)
