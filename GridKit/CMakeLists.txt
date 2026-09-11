configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Definitions.hpp.in
               ${CMAKE_CURRENT_BINARY_DIR}/Definitions.hpp)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/Definitions.hpp DESTINATION include/GridKit)
add_library(definitions INTERFACE)
target_include_directories(
  definitions
  INTERFACE $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}> $<INSTALL_INTERFACE:include>)
add_library(GridKit::definitions ALIAS definitions)
install(TARGETS definitions EXPORT gridkit-targets)

add_library(gridkit_magic_enum INTERFACE)
target_include_directories(
  gridkit_magic_enum
  INTERFACE $<BUILD_INTERFACE:${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include>
            $<INSTALL_INTERFACE:include>)
set_target_properties(
  gridkit_magic_enum
  PROPERTIES EXPORT_NAME magic_enum)
add_library(GridKit::magic_enum ALIAS gridkit_magic_enum)
install(TARGETS gridkit_magic_enum EXPORT gridkit-targets)
install(DIRECTORY ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include/magic_enum DESTINATION include)
install(FILES ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/LICENSE
        DESTINATION share/licenses/GridKit/magic_enum)

# General Utilities and File IO
add_subdirectory(Utilities)

# Memory management utilities
add_subdirectory(MemoryUtilities)

# Create component models
add_subdirectory(Model)

# Local Sparse matrix operations
add_subdirectory(LinearAlgebra)

# Dependency tracking
add_subdirectory(AutomaticDifferentiation)

# Create solvers
add_subdirectory(Solver)

# Testing library
add_subdirectory(Testing)

install(FILES Constants.hpp CommonMath.hpp ScalarTraits.hpp DESTINATION include/GridKit)
