#
#[[

Finds Ipopt include directory and libraries and exports target `IPOPT`

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
mark_as_advanced(IPOPT_LIBRARY IPOPT_INCLUDE_DIR)

if(Ipopt_FOUND AND NOT TARGET IPOPT)
  add_library(IPOPT INTERFACE IMPORTED)
  target_link_libraries(IPOPT INTERFACE ${IPOPT_LIBRARY})
  target_include_directories(IPOPT INTERFACE ${IPOPT_INCLUDE_DIR})
endif()
