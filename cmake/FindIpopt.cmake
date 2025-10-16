# 
#[[

Finds Ipopt include directory and libraries and exports target `Ipopt`

User may set:
- Ipopt_DIR
- IPOPT_ROOT_DIR
- IPOPT_DIR

Author(s):
- Cameron Rutherford <cameron.rutherford@pnnl.gov>

]]

find_library(Ipopt_LIBRARY
  NAMES
  ipopt
  PATHS
  ${IPOPT_DIR} $ENV{IPOPT_DIR} ${IPOPT_ROOT_DIR} ${IPOPT_LIBRARY_DIR}
  ENV LD_LIBRARY_PATH ENV DYLD_LIBRARY_PATH
  PATH_SUFFIXES
  lib64 lib)

if(Ipopt_LIBRARY)
  get_filename_component(IPOPT_LIBRARY_DIR ${Ipopt_LIBRARY} DIRECTORY)
  get_filename_component(IPOPT_DIR ${IPOPT_LIBRARY_DIR} DIRECTORY)
endif()

find_path(Ipopt_INCLUDE_DIR
  NAMES
  IpTNLP.hpp
  PATHS
  ${IPOPT_DIR} ${IPOPT_ROOT_DIR} $ENV{IPOPT_DIR} ${IPOPT_LIBRARY_DIR}/..
  PATH_SUFFIXES
  include
  include/coin
  include/coin-or
  include/coinor)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Ipopt REQUIRED_VARS
    Ipopt_LIBRARY
    Ipopt_INCLUDE_DIR
)

add_library(Ipopt INTERFACE IMPORTED)
target_link_libraries(Ipopt INTERFACE ${IPOPT_LIBRARY})
target_include_directories(Ipopt INTERFACE ${IPOPT_INCLUDE_DIR})

set(Ipopt_DIR ${IPOPT_DIR} CACHE PATH "" FORCE)
mark_as_advanced(Ipopt_INCLUDE_DIR)
mark_as_advanced(Ipopt_LIBRARY)
