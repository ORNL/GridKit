#
#[[

Finds the Ipopt include directory and library and exports target `IPOPT`.

User may set:
- IPOPT_DIR
- IPOPT_ROOT_DIR

Author(s):
- Cameron Rutherford <cameron.rutherford@pnnl.gov>

]]

if(TARGET IPOPT)
  set(Ipopt_FOUND TRUE)
  return()
endif()

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

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Ipopt REQUIRED_VARS IPOPT_LIBRARY IPOPT_INCLUDE_DIR)

if(Ipopt_FOUND AND NOT TARGET IPOPT)
  message(STATUS "Found Ipopt include: ${IPOPT_INCLUDE_DIR}")
  mark_as_advanced(IPOPT_INCLUDE_DIR)
  add_library(IPOPT UNKNOWN IMPORTED)
  set_target_properties(IPOPT PROPERTIES IMPORTED_LOCATION "${IPOPT_LIBRARY}"
                                         INTERFACE_INCLUDE_DIRECTORIES "${IPOPT_INCLUDE_DIR}")
elseif(NOT Ipopt_FOUND)
  if(NOT IPOPT_ROOT_DIR)
    message(STATUS "Ipopt dir not found! Please provide correct filepath.")
    set(IPOPT_DIR
        ${IPOPT_DIR}
        CACHE PATH "Path to Ipopt installation root.")
    unset(IPOPT_INCLUDE_DIR CACHE)
    unset(IPOPT_LIBRARY CACHE)
    unset(IPOPT_LIBRARY_DIR CACHE)
  elseif(NOT IPOPT_LIBRARY)
    message(STATUS "Ipopt library not found! Please provide correct filepath.")
  endif()
  if(IPOPT_ROOT_DIR AND NOT IPOPT_INCLUDE_DIR)
    message(STATUS "Ipopt include directory  not found! Please provide correct path.")
  endif()
endif()
