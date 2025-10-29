# The binary dir is already a global include directory
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/Definitions.hpp.in
  ${CMAKE_BINARY_DIR}/Definitions.hpp)
install(
  FILES ${CMAKE_BINARY_DIR}/Definitions.hpp
  DESTINATION include/GridKit)

# Create component models
add_subdirectory(Model)

# General Utilities and File IO
add_subdirectory(Utilities)

# Local Sparse matrix operations
add_subdirectory(LinearAlgebra)

# Dependency tracking
add_subdirectory(AutomaticDifferentiation)

# Create solvers
add_subdirectory(Solver)

install(
  FILES
    ScalarTraits.hpp
  DESTINATION
    include/GridKit)
