# - Try to find ZMQ
# Once done this will define
# ZMQ_FOUND - System has ZMQ
# ZMQ_INCLUDE_DIRS - The ZMQ include directories
# ZMQ_LIBRARIES - The libraries needed to use ZMQ
# ZMQ_DEFINITIONS - Compiler switches required for using ZMQ

find_path(ZMQ_INCLUDE_DIR zmq.h)
find_library(ZMQ_LIBRARY NAMES zmq)

set(ZMQ_LIBRARIES ${ZMQ_LIBRARY})
set(ZMQ_INCLUDE_DIRS ${ZMQ_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ZMQ REQUIRED_VARS ZMQ_LIBRARIES ZMQ_INCLUDE_DIRS)

add_library(ZMQ UNKNOWN IMPORTED)
set_target_properties(
  ZMQ
  PROPERTIES IMPORTED_LOCATION ${ZMQ_LIBRARIES} INTERFACE_INCLUDE_DIRECTORIES ${ZMQ_INCLUDE_DIRS})
