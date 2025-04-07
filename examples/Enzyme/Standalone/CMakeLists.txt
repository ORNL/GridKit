enzyme_add_executable(
  NAME EnzymeStandaloneScalarCheck
  SOURCES EnzymeScalar.cpp
)

enzyme_add_executable(
  NAME EnzymeStandaloneVectorCheck
  SOURCES EnzymeVector.cpp
  LINK_LIBRARIES GRIDKIT::DenseMatrix
)

enzyme_add_executable(
  NAME EnzymeStandaloneSparseCheck
  SOURCES EnzymeSparse.cpp
  LINK_LIBRARIES GRIDKIT::SparseMatrix
)

add_test(NAME "EnzymeStandaloneScalarCheck" COMMAND ${CMAKE_CURRENT_BINARY_DIR}/EnzymeStandaloneScalarCheck)
add_test(NAME "EnzymeStandaloneVectorCheck" COMMAND ${CMAKE_CURRENT_BINARY_DIR}/EnzymeStandaloneVectorCheck)
add_test(NAME "EnzymeStandaloneSparseCheck" COMMAND ${CMAKE_CURRENT_BINARY_DIR}/EnzymeStandaloneSparseCheck)
