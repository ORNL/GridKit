
# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#]]

#
# gridkit_add_library(<target> SOURCES <file>...
#   [HEADERS <file>...]
#   [LINK_LIBRARIES <item>...]
#   [INCLUDE_DIRECTORIES <dir>...]
#   [COMPILE_OPTIONS <op>...]
#   [OUTPUT_NAME <name>])
#
# Arguments
# - `SOURCES`: source files
# - `HEADERS`: header files that should be installed
# - `LINK_LIBRARIES`: library dependencies
# - `INCLUDE_DIRECTORIES`: header locations
# - `COMPILE_OPTIONS`: extra compiler options
# - `OUTPUT_NAME`: custom library file name (default: `gridkit_<target>`)
#
# loosely based on
# https://github.com/LLNL/sundials/blob/master/cmake/macros/SundialsAddLibrary.cmake
#

macro(gridkit_add_library target)

  set(options "")
  set(oneValueArgs OUTPUT_NAME)
  set(multiValueArgs
      SOURCES
      HEADERS
      LINK_LIBRARIES
      INCLUDE_DIRECTORIES
      COMPILE_OPTIONS)

  # parse arguments
  cmake_parse_arguments(gridkit_add_library
    "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  # add library with sources
  add_library(${target} ${gridkit_add_library_SOURCES})

  # add link libraries
  if(gridkit_add_library_LINK_LIBRARIES)
    target_link_libraries(${target} ${gridkit_add_library_LINK_LIBRARIES})
  endif()

  # set include dirs
  target_include_directories(${target} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
    $<INSTALL_INTERFACE:include/GridKit>)

  if(gridkit_add_library_INCLUDE_DIRECTORIES)
    target_include_directories(${target}
        ${gridkit_add_library_INCLUDE_DIRECTORIES})
  endif()

  # add compile options
  if(gridkit_add_library_COMPILE_OPTIONS)
    target_compile_options(${target} PUBLIC ${gridkit_add_library_COMPILE_OPTIONS})
  endif()

  # add alias library
  add_library(GRIDKIT::${target} ALIAS ${target})

  # set output name
  if(gridkit_add_library_OUTPUT_NAME)
    set(_output_name ${gridkit_add_library_OUTPUT_NAME})
  else()
    set(_output_name gridkit_${target})
  endif()
  set_target_properties(${target} PROPERTIES OUTPUT_NAME ${_output_name})

  # set the library version
  set_target_properties(${target} PROPERTIES
    VERSION ${PACKAGE_VERSION}
    SOVERSION ${PACKAGE_VERSION_MAJOR})

  # install
  install(TARGETS ${target} DESTINATION lib EXPORT gridkit-targets)
  if(gridkit_add_library_HEADERS)
    cmake_path(RELATIVE_PATH CMAKE_CURRENT_SOURCE_DIR
      BASE_DIRECTORY ${CMAKE_SOURCE_DIR}/src
      OUTPUT_VARIABLE _rel_path)
    install(FILES ${gridkit_add_library_HEADERS}
      DESTINATION include/GridKit/${_rel_path})
  endif()

endmacro()
