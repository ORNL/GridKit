#
# Copyright (c) 2017, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by Slaven Peles <peles2@llnl.gov>.
# LLNL-CODE-718378.
# All rights reserved.
#
# This file is part of GridKit™. For details, see github.com/LLNL/GridKit
# Please also read the LICENSE file.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# - Redistributions of source code must retain the above copyright notice,
#   this list of conditions and the disclaimer below.
# - Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the disclaimer (as noted below) in the
#   documentation and/or other materials provided with the distribution.
# - Neither the name of the LLNS/LLNL nor the names of its contributors may
#   be used to endorse or promote products derived from this software without
#   specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
# CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
# INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL
# SECURITY, LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
# OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISINGIN ANY
# WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
# THE POSSIBILITY OF SUCH DAMAGE.
#
# Lawrence Livermore National Laboratory is operated by Lawrence Livermore
# National Security, LLC, for the U.S. Department of Energy, National
# Nuclear Security Administration under Contract DE-AC52-07NA27344.
#
# This document was prepared as an account of work sponsored by an agency
# of the United States government. Neither the United States government nor
# Lawrence Livermore National Security, LLC, nor any of their employees
# makes any warranty, expressed or implied, or assumes any legal liability
# or responsibility for the accuracy, completeness, or usefulness of any
# information, apparatus, product, or process disclosed, or represents that
# its use would not infringe privately owned rights. Reference herein to
# any specific commercial product, process, or service by trade name,
# trademark, manufacturer, or otherwise does not necessarily constitute or
# imply its endorsement, recommendation, or favoring by the United States
# government or Lawrence Livermore National Security, LLC. The views and
# opinions of authors expressed herein do not necessarily state or reflect
# those of the United States government or Lawrence Livermore National
# Security, LLC, and shall not be used for advertising or product
# endorsement purposes.
#
#
# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#    - Slaven Peles <peless@ornl.gov>
# ]]

cmake_minimum_required(VERSION 3.12)

project(gridkit)

set(PACKAGE_NAME "GRIDKIT")
set(PACKAGE_STRING "GRIDKIT 0.0.7")
set(PACKAGE_TARNAME "gridkit")

set(PACKAGE_VERSION_MAJOR "0")
set(PACKAGE_VERSION_MINOR "0")
set(PACKAGE_VERSION_PATCH "7")

set(PACKAGE_VERSION "${PACKAGE_VERSION_MAJOR}.${PACKAGE_VERSION_MINOR}.${PACKAGE_VERSION_PATCH}")

# Ipopt support is disabled by default
option(GRIDKIT_ENABLE_IPOPT "Enable Ipopt support" ON)

# SUNDIALS support is enabled by default
option(GRIDKIT_ENABLE_SUNDIALS "Enable SUNDIALS support" ON)

# Enable KLU
option(GRIDKIT_ENABLE_SUNDIALS_SPARSE "Enable SUNDIALS sparse linear solvers" ON)


set(CMAKE_MACOSX_RPATH 1)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMake)

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


# TODO: Probably beter to set a debug interface target
set(CMAKE_CXX_FLAGS_DEBUG "-Wall -O0 -g -DDEBUG")

set(CMAKE_CXX_STANDARD 17)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

option(GRIDKIT_BUILD_SHARED "Build shared libraries" ON)
option(GRIDKIT_BUILD_STATIC "Build static libraries" OFF)

if(GRIDKIT_ENABLE_IPOPT)
  include(FindIpopt)
endif()
if(GRIDKIT_ENABLE_SUNDIALS)
  find_package(SUNDIALS 6.0.0 REQUIRED CONFIG 
               PATHS ${SUNDIALS_DIR}
                     ${SUNDIALS_DIR}/lib/cmake/sundials)
  message(STATUS "SUNDIALS configuration found: ${SUNDIALS_CONFIG}")
endif()
if(GRIDKIT_ENABLE_SUNDIALS_SPARSE)
  include(FindSuiteSparse)
endif()

# Macro that adds libraries
include(GridkitAddLibrary)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Create component models
add_subdirectory(ComponentLib)

# General Utilities and File IO
add_subdirectory(Utilities)

#Local Sparse matrix operations
add_subdirectory(SparseMatrix)

# Create solvers
add_subdirectory(Solver)

# Create examples and tests
enable_testing()
add_subdirectory(Examples)
  
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
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/CMake/Config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/GridKitConfig.cmake
  INSTALL_DESTINATION share/cmake/gridkit
  NO_SET_AND_CHECK_MACRO
  NO_CHECK_REQUIRED_COMPONENTS_MACRO)

# Install configuration file
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/GridKitConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/GridKitConfigVersion.cmake
  DESTINATION share/cmake/gridkit)
