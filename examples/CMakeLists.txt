# [[
#  Author(s):
#    - Slaven Peles <peless@ornl.gov>
#]]

set(GRIDKIT_EXAMPLES_INSTALL_ROOT "share/gridkit/examples")

#
# gridkit_example_current_install_path(<output-variable>)
#
# Get the full installation destination for the current example directory.
#
# Result is stored in <output-variable>
#
macro(gridkit_example_current_install_path _out)
  cmake_path(
    RELATIVE_PATH
    CMAKE_CURRENT_SOURCE_DIR
    BASE_DIRECTORY
    ${CMAKE_SOURCE_DIR}/examples
    OUTPUT_VARIABLE
    _rel_path)
  set(${_out} ${GRIDKIT_EXAMPLES_INSTALL_ROOT}/${_rel_path})
endmacro()

#
# gridkit_example_add_file(<filename>)
#
# This macro will add a file (e.g., JSON) that is needed for an example to the
# build tree and the install tree following the same directory structure as the
# source tree under `examples`.
#
# <filename> may be an absolute path or a path relative to
# CMAKE_CURRENT_SOURCE_DIR.
#
macro(gridkit_example_add_file _filename)
  get_filename_component(__fname ${_filename} NAME)
  get_filename_component(
    _source_filepath
    ${_filename}
    ABSOLUTE
    BASE_DIR
    ${CMAKE_CURRENT_SOURCE_DIR})
  if(NOT EXISTS ${_source_filepath})
    message(FATAL_ERROR "File does not exist: ${_source_filepath}")
  endif()
  set(_binary_filepath ${CMAKE_CURRENT_BINARY_DIR}/${__fname})
  configure_file(${_source_filepath} ${_binary_filepath} COPYONLY)
  gridkit_example_current_install_path(_ex_install_path)
  install(FILES ${_source_filepath} DESTINATION ${_ex_install_path})
endmacro()

add_subdirectory(PowerFlow)
add_subdirectory(PowerElectronics)
add_subdirectory(Experimental)

if(TARGET SUNDIALS::idas)
  add_subdirectory(PhasorDynamics)
  add_subdirectory(Network)
endif()

add_subdirectory(Consumer)
