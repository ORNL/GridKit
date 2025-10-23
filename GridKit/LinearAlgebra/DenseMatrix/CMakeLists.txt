
add_library(DenseMatrix INTERFACE)
target_include_directories(DenseMatrix INTERFACE ${CMAKE_SOURCE_DIR})

add_library(GRIDKIT::DenseMatrix ALIAS DenseMatrix)
