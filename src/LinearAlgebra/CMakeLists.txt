
add_subdirectory(SparseMatrix)
add_subdirectory(DenseMatrix)

gridkit_add_library(linear_workspace
    SOURCES
      LinAlgWorkspaceCpu.cpp)
