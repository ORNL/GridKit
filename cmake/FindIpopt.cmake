#
#[[

Finds Ipopt include directory and libraries and exports target `Ipopt`

User may set:
- IPOPT_ROOT_DIR

Author(s):
- Cameron Rutherford <cameron.rutherford@pnnl.gov>

]]

find_library(
  IPOPT_LIBRARY
  NAMES ipopt
  PATHS ${IPOPT_DIR}
        $ENV{IPOPT_DIR}
        ${IPOPT_ROOT_DIR}
        ${IPOPT_LIBRARY_DIR}
        ENV
        LD_LIBRARY_PATH
        ENV
        DYLD_LIBRARY_PATH
  PATH_SUFFIXES lib64 lib)

if(IPOPT_LIBRARY)
  set(IPOPT_LIBRARY CACHE FILEPATH "Path to Ipopt library")
  message(STATUS "Found Ipopt library: " ${IPOPT_LIBRARY})
  get_filename_component(
    IPOPT_LIBRARY_DIR
    ${IPOPT_LIBRARY}
    DIRECTORY
    CACHE
    "Ipopt library directory")
  mark_as_advanced(IPOPT_LIBRARY IPOPT_LIBRARY_DIR)
  if(NOT IPOPT_DIR)
    get_filename_component(
      IPOPT_DIR
      ${IPOPT_LIBRARY_DIR}
      DIRECTORY
      CACHE)
  endif()
endif()

find_path(
  IPOPT_INCLUDE_DIR
  NAMES IpTNLP.hpp
  PATHS ${IPOPT_DIR}
        ${IPOPT_ROOT_DIR}
        $ENV{IPOPT_DIR}
        ${IPOPT_LIBRARY_DIR}/..
  PATH_SUFFIXES
    include
    include/coin
    include/coin-or
    include/coinor)

if(IPOPT_LIBRARY)
  message(STATUS "Found Ipopt include: ${IPOPT_INCLUDE_DIR}")
  mark_as_advanced(IPOPT_INCLUDE_DIR)
  add_library(IPOPT INTERFACE IMPORTED)
  target_link_libraries(IPOPT INTERFACE ${IPOPT_LIBRARY})
  target_include_directories(IPOPT INTERFACE ${IPOPT_INCLUDE_DIR})
else()
  if(NOT IPOPT_ROOT_DIR)
    message(STATUS "Ipopt dir not found! Please provide correct filepath.")
    set(IPOPT_DIR
        ${IPOPT_DIR}
        CACHE PATH "Path to Ipopt installation root.")
    unset(IPOPT_INCLUDE_DIR CACHE)
    unset(IPOPT_LIBRARY CACHE)
    unset(IPOPT_LIBRARY_DIR CACHE)
  elseif(NOT IPOPT_LIB)
    message(STATUS "Ipopt library not found! Please provide correct filepath.")
  endif()
  if(IPOPT_ROOT_DIR AND NOT IPOPT_INCLUDE_DIR)
    message(STATUS "Ipopt include directory  not found! Please provide correct path.")
  endif()
endif()
