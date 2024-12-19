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

#[[

Finds Enzyme Clang plugin

User may set:
- ENZYME_DIR

Author(s):
- Asher Mancinelli <ashermancinelli@gmail.com>

]]

if(NOT ${CMAKE_C_COMPILER_ID} STREQUAL "Clang" OR NOT ${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
  message(FATAL_ERROR "Enzyme may only be enabled if building with Clang")
endif()

find_library(ENZYME_LLVM_PLUGIN_LIBRARY
  NAMES
  LLVMEnzyme-15.so
  LLVMEnzyme-14.so
  LLVMEnzyme-13.so
  LLVMEnzyme-12.so
  LLVMEnzyme-11.so
  LLVMEnzyme-10.so
  LLVMEnzyme-9.so
  PATHS
  ${ENZYME_DIR} $ENV{ENZYME_DIR}
  ENV LD_LIBRARY_PATH ENV DYLD_LIBRARY_PATH
  PATH_SUFFIXES
  lib64 lib)

find_library(ENZYME_CLANG_PLUGIN_LIBRARY
  NAMES
  ClangEnzyme-15.so
  ClangEnzyme-14.so
  ClangEnzyme-13.so
  ClangEnzyme-12.so
  ClangEnzyme-11.so
  ClangEnzyme-10.so
  ClangEnzyme-9.so
  PATHS
  ${ENZYME_DIR} $ENV{ENZYME_DIR}
  ENV LD_LIBRARY_PATH ENV DYLD_LIBRARY_PATH
  PATH_SUFFIXES
  lib64 lib)

find_program(GRIDKIT_LLVM_LINK llvm-link)
find_program(GRIDKIT_OPT opt)

if("${GRIDKIT_LLVM_LINK}" STREQUAL "GRIDKIT_LLVM_LINK-NOTFOUND")
  message(FATAL_ERROR "Could not find llvm-link! Will not be able to build Enzyme targets.")
endif()

if("${GRIDKIT_OPT}" STREQUAL "GRIDKIT_OPT-NOTFOUND")
  message(FATAL_ERROR "Could not find opt! Will not be able to build Enzyme targets.")
endif()

message(STATUS "${ENZYME_CLANG_PLUGIN_LIBRARY};${ENZYME_LLVM_PLUGIN_LIBRARY}")
if((NOT ${ENZYME_LLVM_PLUGIN_LIBRARY}) OR (NOT ${ENZYME_CLANG_PLUGIN_LIBRARY}))
  set(ENZYME_CLANG_PLUGIN_LIBRARY CACHE FILEPATH "Path to Enzyme Clang plugin library")
  set(ENZYME_LLVM_PLUGIN_LIBRARY CACHE FILEPATH "Path to Enzyme LLVM plugin library")
  message(STATUS "Found Enzyme clang plugin: ${ENZYME_CLANG_PLUGIN_LIBRARY}")
  message(STATUS "Found Enzyme LLVM plugin: ${ENZYME_LLVM_PLUGIN_LIBRARY}")
  get_filename_component(ENZYME_LIBRARY_DIR ${ENZYME_CLANG_PLUGIN_LIBRARY} DIRECTORY CACHE "Enzyme library directory")
else()
  message(FATAL_ERROR "Enzyme library not found!"
    " Set ENZYME_DIR to Enzyme's installation prefix.")
endif()

macro(enzyme_add_executable)
  set(options)
  set(oneValueArgs NAME)
  set(multiValueArgs SOURCES LINK_LIBRARIES)
  cmake_parse_arguments(enzyme_add_executable "${options}" "${oneValueArgs}"
    "${multiValueArgs}" ${ARGN})

  set(PHASE2 "${CMAKE_CURRENT_BINARY_DIR}/${enzyme_add_executable_NAME}.bc")
  set(PHASE3 "${CMAKE_CURRENT_BINARY_DIR}/${enzyme_add_executable_NAME}_enzyme.ll")
  set(PHASE4 "${CMAKE_CURRENT_BINARY_DIR}/${enzyme_add_executable_NAME}_opt.ll")
  set(PHASE5 "${CMAKE_CURRENT_BINARY_DIR}/${enzyme_add_executable_NAME}")

  set(OBJS "")

  foreach(SRC ${enzyme_add_executable_SOURCES})
    set(PHASE0 "${CMAKE_CURRENT_SOURCE_DIR}/${SRC}")
    set(PHASE1 "${CMAKE_CURRENT_BINARY_DIR}/${enzyme_add_executable_NAME}_${SRC}_compile.o")
    add_custom_command(
      DEPENDS ${PHASE0}
      OUTPUT ${PHASE1}
      COMMAND ${CMAKE_CXX_COMPILER} -flto -c ${PHASE0} -O2 -fno-vectorize -ffast-math -fno-unroll-loops -o ${PHASE1}
      COMMENT "Compiling ${SRC} to object file for target ${enzyme_add_executable_NAME}"
      )
    set(OBJS "${OBJS} ${PHASE1}")
  endforeach()

  cmake_language(EVAL CODE "
  add_custom_command(
    DEPENDS ${OBJS}
    OUTPUT ${PHASE2}
    COMMAND ${GRIDKIT_LLVM_LINK} ${OBJS} -o ${PHASE2}
    COMMENT \"Linking object files to LLVM bytecode for target ${enzyme_add_executable_NAME}\"
    )
  ")

  add_custom_command(
    DEPENDS ${PHASE2}
    OUTPUT ${PHASE3}
    COMMAND ${GRIDKIT_OPT} ${PHASE2} --enable-new-pm=0 -load=${ENZYME_LLVM_PLUGIN_LIBRARY} -enzyme -o ${PHASE3} -S
    COMMENT "Running Enzyme opt pass on target ${enzyme_add_executable_NAME}"
    )

  add_custom_command(
    DEPENDS ${PHASE3}
    OUTPUT ${PHASE4}
    COMMAND ${GRIDKIT_OPT} ${PHASE3} -O2 -o ${PHASE4} -S
    COMMENT "Running remaining opt passes on target ${enzyme_add_executable_NAME}"
    )

  add_custom_command(
    DEPENDS ${PHASE4}
    OUTPUT ${PHASE5}
    COMMAND ${CMAKE_CXX_COMPILER} ${PHASE4} -o ${PHASE5}
    )

  add_custom_target(
    ${enzyme_add_executable_NAME} ALL
    DEPENDS ${PHASE5}
    )
endmacro()
