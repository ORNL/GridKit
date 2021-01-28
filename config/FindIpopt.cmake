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

Finds Ipopt include directory and libraries and exports target `Ipopt`

User may set:
- IPOPT_ROOT_DIR

]]

find_library(IPOPT_LIBRARY
  NAMES
  ipopt
  PATHS
  ${IPOPT_DIR} $ENV{IPOPT_DIR} ${IPOPT_ROOT_DIR}
  ENV LD_LIBRARY_PATH ENV DYLD_LIBRARY_PATH
  PATH_SUFFIXES
  lib64 lib)

if(IPOPT_LIBRARY)
  set(IPOPT_LIBRARY CACHE FILEPATH "Path to Ipopt library")
  message(STATUS "Found Ipopt library: " ${IPOPT_LIBRARY})
  get_filename_component(IPOPT_LIBRARY_DIR ${IPOPT_LIBRARY} DIRECTORY CACHE)
  set(IPOPT_LIBRARY_DIR CACHE PATH "Path to Ipopt library")
  if(NOT IPOPT_DIR)
    get_filename_component(IPOPT_DIR ${IPOPT_LIBRARY_DIR} DIRECTORY CACHE)
  endif()
endif()

find_path(IPOPT_INCLUDE_DIR
  NAMES
  IpTNLP.hpp
  PATHS
  ${IPOPT_DIR} ${IPOPT_ROOT_DIR} $ENV{IPOPT_DIR} ${IPOPT_LIBRARY_DIR}/..
  PATH_SUFFIXES
  include
  include/coin
  include/coin-or
  include/coinor)

if(IPOPT_LIBRARY AND IPOPT_INCLUDE_DIR)
  set(IPOPT_INCLUDE_DIR CACHE PATH "Path to Ipopt header files")
  message(STATUS "Found Ipopt include directory: " ${IPOPT_INCLUDE_DIR})
  add_library(Ipopt INTERFACE)
  target_link_libraries(Ipopt INTERFACE ${IPOPT_LIBRARY})
  target_include_directories(Ipopt INTERFACE ${IPOPT_INCLUDE_DIR})
else()
  if(NOT IPOPT_ROOT_DIR)
    message(STATUS "Ipopt dir not found! Please provide correct filepath.")
    set(IPOPT_DIR CACHE PATH "Path to Ipopt installation root.")
    unset(IPOPT_INCLUDE_DIR CACHE)
    unset(IPOPT_LIBRARY CACHE)
  elseif(NOT IPOPT_LIB)
    message(STATUS "Ipopt library not found! Please provide correct filepath.")
  endif()
  if(IPOPT_ROOT_DIR AND NOT IPOPT_INCLUDE_DIR)
    message(STATUS "Ipopt include directory  not found! Please provide correct path.")
  endif()
endif()

