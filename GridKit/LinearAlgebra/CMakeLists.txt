add_subdirectory(SparseMatrix)
add_subdirectory(Vector)
add_subdirectory(Solver)

install(FILES DenseExpression.hpp SparseJet.hpp DESTINATION include/GridKit/LinearAlgebra)
