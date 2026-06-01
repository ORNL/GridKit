add_subdirectory(SparseMatrix)
add_subdirectory(DenseMatrix)

install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/MemoryUtils.hpp DESTINATION include/GridKit/LinearAlgebra)
