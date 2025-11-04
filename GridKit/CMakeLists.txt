configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/Definitions.hpp.in
  ${CMAKE_CURRENT_BINARY_DIR}/Definitions.hpp)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/Definitions.hpp
  DESTINATION include/GridKit)
add_library(definitions INTERFACE)
target_include_directories(definitions INTERFACE
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
  $<INSTALL_INTERFACE:include>)
add_library(GridKit::definitions ALIAS definitions)
install(TARGETS definitions EXPORT gridkit-targets)

# General Utilities and File IO
add_subdirectory(Utilities)

# Create component models
add_subdirectory(Model)

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
