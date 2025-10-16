#
#
# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#    - Slaven Peles <peless@ornl.gov>
# ]]

cmake_minimum_required(VERSION 3.13)

project(gridkit)

set(PACKAGE_NAME "GRIDKIT")
set(PACKAGE_STRING "GRIDKIT 0.1.0")
set(PACKAGE_TARNAME "gridkit")

set(PACKAGE_VERSION_MAJOR "0")
set(PACKAGE_VERSION_MINOR "1")
set(PACKAGE_VERSION_PATCH "0")

set(GRIDKIT_THIRD_PARTY_DIR ${PROJECT_SOURCE_DIR}/third-party)

set(PACKAGE_VERSION "${PACKAGE_VERSION_MAJOR}.${PACKAGE_VERSION_MINOR}.${PACKAGE_VERSION_PATCH}")

set(CMAKE_CXX_STANDARD 20)

if (MSVC)
  set(CMAKE_CXX_FLAGS "/Wall")
else()
  set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wconversion -Wpedantic")
endif()

# Ipopt support is disabled by default
option(GRIDKIT_ENABLE_IPOPT "Enable Ipopt support" OFF)

# SUNDIALS support is disabled by default
option(GRIDKIT_ENABLE_SUNDIALS "Enable SUNDIALS support" OFF)

# Enzyme support is disabled by default
option(GRIDKIT_ENABLE_ENZYME "Enable automatic differentiation with Enzyme" OFF)

# Enable the address sanitizer
option(GRIDKIT_ENABLE_ASAN "Enable the address sanitizer" OFF)

if(GRIDKIT_ENABLE_ASAN)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address")
endif()

# Enable the undefined behavior sanitizer
option(GRIDKIT_ENABLE_UBSAN "Enable the undefined behavior sanitizer" OFF)

if(GRIDKIT_ENABLE_UBSAN)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=undefined -fno-omit-frame-pointer")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=undefined")
endif()

set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS 1)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# https://gitlab.kitware.com/cmake/community/-/wikis/doc/cmake/RPATH-handling#always-full-rpath
# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used when installing, but only if it's not a system directory
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif("${isSystemDir}" STREQUAL "-1")


include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src ${CMAKE_BINARY_DIR})

option(GRIDKIT_BUILD_SHARED "Build shared libraries" ON)
option(GRIDKIT_BUILD_STATIC "Build static libraries" OFF)

if(GRIDKIT_ENABLE_IPOPT)
  find_package(Ipopt REQUIRED)
  enable_language(Fortran) # Needed for linking to HSL
endif()
if(GRIDKIT_ENABLE_SUNDIALS)
  list(APPEND CMAKE_PREFIX_PATH ${SUNDIALS_DIR})
  find_package(SUNDIALS 7.4.0 REQUIRED)
  message(STATUS "SUNDIALS configuration found: ${SUNDIALS_CONFIG}")
  if(TARGET SUNDIALS::sunlinsolklu)
    set(GRIDKIT_ENABLE_SUNDIALS_SPARSE ON CACHE INTERNAL BOOL "Sparse solver enabled in SUNDIALS.")
    message(STATUS "KLU support enabled in SUNDIALS.")
  endif()
endif()

if(GRIDKIT_ENABLE_ENZYME)
  include(FindEnzyme)
endif()

# The binary dir is already a global include directory
configure_file(
  ${CMAKE_SOURCE_DIR}/src/Definitions.hpp.in
  ${CMAKE_BINARY_DIR}/Definitions.hpp)
install(
  FILES ${CMAKE_BINARY_DIR}/Definitions.hpp
  DESTINATION include
  )

# Macro that adds libraries
include(GridkitAddLibrary)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

add_subdirectory(src)

# Create examples and tests
enable_testing()
add_subdirectory(examples)
add_subdirectory(tests)

export(EXPORT gridkit-targets FILE ${CMAKE_CURRENT_BINARY_DIR}/GridKitTargets.cmake)

# Configuring exporting cmake config
install(EXPORT gridkit-targets
  FILE GridKitTargets.cmake
  NAMESPACE GRIDKIT::
  DESTINATION share/cmake/gridkit)

include(CMakePackageConfigHelpers)

# Basic version file
write_basic_package_version_file(
  GridKitConfigVersion.cmake
  VERSION ${PACKAGE_VERSION}
  COMPATIBILITY SameMajorVersion)

# Generate config file that includes exports
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/Config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/GridKitConfig.cmake
  INSTALL_DESTINATION share/cmake/gridkit
  NO_SET_AND_CHECK_MACRO
  NO_CHECK_REQUIRED_COMPONENTS_MACRO)

# Install configuration file
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/GridKitConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/GridKitConfigVersion.cmake
  DESTINATION share/cmake/gridkit)
