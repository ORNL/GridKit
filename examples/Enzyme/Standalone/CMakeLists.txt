add_executable(EnzymeStandaloneScalarCheck EnzymeScalar.cpp)
target_link_libraries(
  EnzymeStandaloneScalarCheck ClangEnzymeFlags GridKit::Utilities)

add_executable(EnzymeStandaloneVectorCheck EnzymeVector.cpp)
target_link_libraries(
  EnzymeStandaloneVectorCheck
  ClangEnzymeFlags
  GridKit::DenseMatrix
  GridKit::Utilities)

add_executable(EnzymeStandaloneSparseCheck EnzymeSparse.cpp)
target_compile_options(EnzymeStandaloneSparseCheck PUBLIC -mllvm -enzyme-auto-sparsity=1)
target_link_libraries(
  EnzymeStandaloneSparseCheck
  ClangEnzymeFlags
  GridKit::sparse_matrix
  GridKit::Utilities)

add_test(NAME "EnzymeStandaloneScalarCheck"
         COMMAND ${CMAKE_CURRENT_BINARY_DIR}/EnzymeStandaloneScalarCheck)
add_test(NAME "EnzymeStandaloneVectorCheck"
         COMMAND ${CMAKE_CURRENT_BINARY_DIR}/EnzymeStandaloneVectorCheck)
add_test(NAME "EnzymeStandaloneSparseCheck"
         COMMAND ${CMAKE_CURRENT_BINARY_DIR}/EnzymeStandaloneSparseCheck)
