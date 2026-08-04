gridkit_add_library(
  optimization_core
  INTERFACE_TARGET
  HEADERS DerivativeStructure.hpp Evaluator.hpp SparseMatrixAssembler.hpp
  LINK_LIBRARIES
    INTERFACE
    GridKit::dense_vector
    INTERFACE
    GridKit::sparse_matrix)

if(GRIDKIT_ENABLE_ENZYME)
  gridkit_add_library(
    optimization_enzyme
    INTERFACE_TARGET
    LINK_LIBRARIES
      INTERFACE
      GridKit::optimization_core
      INTERFACE
      ClangEnzymeFlags
    COMPILE_OPTIONS INTERFACE -fno-math-errno)
endif()
