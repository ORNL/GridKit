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

Finds Sutiesparse include directory and libraries and exports target `Suitesparse`

User may set:
- SUITESPARSE_ROOT_DIR

Author(s):
- Cameron Rutherford <cameron.rutherford@pnnl.gov>

]]
set(SUITESPARSE_MODULES
  amd
  colamd
  klu
  suitesparseconfig)

find_library(SUITESPARSE_LIBRARY
  NAMES
  ${SUITESPARSE_MODULES}
  PATHS
  ${SUITESPARSE_DIR} $ENV{SUITESPARSE_DIR} ${SUITESPARSE_ROOT_DIR}
  ENV LD_LIBRARY_PATH ENV DYLD_LIBRARY_PATH
  PATH_SUFFIXES
  lib64 lib)

if(SUITESPARSE_LIBRARY)
  set(SUITESPARSE_LIBRARY CACHE FILEPATH "Path to Suitesparse library")
  get_filename_component(SUITESPARSE_LIBRARY_DIR ${SUITESPARSE_LIBRARY} DIRECTORY CACHE "Suitesparse library directory")
  message(STATUS "Found Suitesparse libraries in: " ${SUITESPARSE_LIBRARY_DIR})
  mark_as_advanced(SUITESPARSE_LIBRARY SUITESPARSE_LIBRARY_DIR)
  if(NOT SUITESPARSE_DIR) 
    get_filename_component(SUITESPARSE_DIR ${SUITESPARSE_LIBRARY_DIR} DIRECTORY CACHE)
  endif()
endif()

# Find SUITESPARSE header path and ensure all needed files are there
find_path(SUITESPARSE_INCLUDE_DIR
  NAMES
  amd.h
  colamd.h
  klu.h
  SuiteSparse_config.h
  PATHS
  ${SUITESPARSE_DIR} $ENV{SUITESPARSE_DIR} ${SUITESPARSE_ROOT_DIR} ${SUITESPARSE_LIBRARY_DIR}/..
  PATH_SUFFIXES
  include)

if(SUITESPARSE_LIBRARY)
  message(STATUS "Found Suitesparse include: ${SUITESPARSE_INCLUDE_DIR}")
  mark_as_advanced(SUITESPARSE_INCLUDE_DIR)
  unset(SUITESPARSE_LIBRARY)
  add_library(SUITESPARSE INTERFACE IMPORTED)
  target_include_directories(SUITESPARSE INTERFACE ${SUITESPARSE_INCLUDE_DIR})
  foreach(mod ${SUITESPARSE_MODULES})
    find_library(suitesparse_${mod} 
      NAMES ${mod}
      HINTS ${SUITESPARSE_LIBRARY_DIR})
    if(suitesparse_${mod})
      message(STATUS "Found suitesparse internal library " ${mod})
      target_link_libraries(SUITESPARSE INTERFACE ${suitesparse_${mod}})
      mark_as_advanced(suitesparse_${mod})
    else()
      message(SEND_ERROR "Suitesparse internal library " ${mod} " not found")
    endif()
  endforeach(mod)
else()
  if(NOT SUITESPARSE_ROOT_DIR)
    message(STATUS "Suitesparse dir not found! Please provide correct filepath.")
    set(SUITESPARSE_DIR ${SUITESPARSE_DIR} CACHE PATH "Path to Suitesparse installation root.")
    unset(SUITESPARSE_LIBRARY CACHE)
    unset(SUITESPARSE_INCLUDE_DIR CACHE)
    unset(SUITESPARSE_LIBRARY_DIR CACHE)
  elseif(NOT SUITESPARSE_LIBRARY)
    message(STATUS "Suitesparse library not found! Please provide correct filepath.")
  endif()
  if(SUITESPARSE_ROOT_DIR AND NOT SUITESPARSE_INCLUDE_DIR)
    message(STATUS "Suitesparse include dir not found! Please provide correct filepath.")
  endif()
endif()
