#
#
# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#    - Slaven Peles <peless@ornl.gov>
# ]]

cmake_minimum_required(VERSION 3.13)

project(GridKit 
    VERSION 0.1.0
    LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if (MSVC)
  set(CMAKE_CXX_FLAGS "/Wall")
else()
  set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wconversion -Wpedantic")
endif()

set(GRIDKIT_THIRD_PARTY_DIR ${PROJECT_SOURCE_DIR}/third-party)

# Ipopt support is disabled by default
option(GridKit_ENABLE_IPOPT "Enable Ipopt support" OFF)

# SUNDIALS support is disabled by default
option(GridKit_ENABLE_SUNDIALS "Enable SUNDIALS support" ON)

# Enzyme support is disabled by default
option(GridKit_ENABLE_ENZYME "Enable automatic differentiation with Enzyme" OFF)

# ReSolve support is disabled by default
option(GridKit_ENABLE_RESOLVE "Enable ReSolve support" OFF)

# Enable the address sanitizer
option(GridKit_ENABLE_ASAN "Enable the address sanitizer" OFF)

# Enable the undefined behavior sanitizer
option(GridKit_ENABLE_UBSAN "Enable the undefined behavior sanitizer" OFF)

# This allows use of "GRIDKIT_*" versions of the above options
list(APPEND _gridkit_enable_options IPOPT SUNDIALS ENZYME RESOLVE ASAN UBSAN)
foreach(_opt IN LISTS _gridkit_enable_options)
  if(NOT DEFINED GRIDKIT_ENABLE_${_opt})
    set(GRIDKIT_ENABLE_${_opt} ${GridKit_ENABLE_${_opt}})
  endif()
endforeach()

# Sanitizers
if(GRIDKIT_ENABLE_ASAN)
  set(GRIDKIT_COMPILE_OPTIONS ${GRIDKIT_COMPILE_OPTIONS} -fsanitize=address -fno-omit-frame-pointer CACHE INTERNAL STRING)
  set(GRIDKIT_LINK_OPTIONS ${GRIDKIT_LINK_OPTIONS} -fsanitize=address CACHE INTERNAL STRING)
endif()

if(GRIDKIT_ENABLE_UBSAN)
  set(GRIDKIT_COMPILE_OPTIONS ${GRIDKIT_COMPILE_OPTIONS} -fsanitize=undefined -fno-omit-frame-pointer CACHE INTERNAL STRING)
  set(GRIDKIT_LINK_OPTIONS ${GRIDKIT_LINK_OPTIONS} -fsanitize=undefined CACHE INTERNAL STRING)
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

option(BUILD_SHARED_LIBS "Build shared libraries" ON)

if(GRIDKIT_ENABLE_IPOPT)
  include(FindIpopt)
  enable_language(Fortran) # Needed for linking to HSL
endif()
if(GRIDKIT_ENABLE_SUNDIALS)
  enable_language(C)
  find_package(SUNDIALS 7.4.0 REQUIRED CONFIG
               PATHS ${SUNDIALS_DIR}
                     ${SUNDIALS_DIR}/lib/cmake/sundials)
  message(STATUS "SUNDIALS configuration found: ${SUNDIALS_CONFIG}")
  if(TARGET SUNDIALS::sunlinsolklu)
    set(GRIDKIT_ENABLE_SUNDIALS_SPARSE ON CACHE INTERNAL BOOL "Sparse solver enabled in SUNDIALS.")
    message(STATUS "KLU support enabled in SUNDIALS.")
  endif()
endif()

if(GRIDKIT_ENABLE_ENZYME)
  include(FindEnzyme)
endif()

if(GRIDKIT_ENABLE_RESOLVE)
  find_package(ReSolve CONFIG 
    PATHS ${ReSolve_DIR} ${ReSolve_DIR}/share/resolve/cmake
    ENV LD_LIBRARY_PATH ENV DYLD_LIBRARY_PATH
    REQUIRED)
endif()

# Macro that adds libraries
include(GridkitAddLibrary)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

add_subdirectory(GridKit)

# Create examples and tests
include(CTest)
add_subdirectory(examples)
add_subdirectory(tests)

# Documentation
add_subdirectory(docs)

export(EXPORT gridkit-targets FILE ${CMAKE_CURRENT_BINARY_DIR}/GridKitTargets.cmake)

# Configuring exporting cmake config
install(EXPORT gridkit-targets
  FILE GridKitTargets.cmake
  NAMESPACE GridKit::
  DESTINATION share/cmake/gridkit)

include(CMakePackageConfigHelpers)

# Basic version file
write_basic_package_version_file(
  GridKitConfigVersion.cmake
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
