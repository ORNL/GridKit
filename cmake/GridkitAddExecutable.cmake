# add_executable macro loosely based on add_library

macro(gridkit_add_executable target)

  set(options TEST)
  set(multiValueArgs SOURCES LINK_LIBRARIES COMPILE_OPTIONS)

  # parse arguments
  cmake_parse_arguments(gridkit_add_executable
    "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  # Before going any further and configuring this library, make sure we can actually link its requirements.
  foreach(_add_lib ${gridkit_add_executable_LINK_LIBRARIES})
    if(NOT (${_add_lib} MATCHES "PUBLIC" OR ${_add_lib} MATCHES "PRIVATE" OR ${_add_lib} MATCHES "INTERFACE"))
      if(NOT TARGET ${_add_lib})
        get_property(_skipped_libs GLOBAL PROPERTY GRIDKIT_SKIPPED_OPTIONAL_LIBS)

        if(${_add_lib} IN_LIST _skipped_libs)
          message(STATUS "Disabling GridKit executable ${target} due to disabled optional library: ${_add_lib}")
          return()
        else()
          message(SEND_ERROR "Failed to configure GridKit executable due to missing library: ${_add_lib}")
        endif()
      endif()
    endif()
  endforeach()

  # -- Create executable --

  add_executable(${target} ${gridkit_add_executable_SOURCES})

  if(gridkit_add_executable_COMPILE_OPTIONS)
    target_compile_options(${target} PUBLIC ${gridkit_add_executable_COMPILE_OPTIONS})
  endif()
  
  if(gridkit_add_executable_LINK_LIBRARIES)
    target_link_libraries(${target} ${gridkit_add_executable_LINK_LIBRARIES})
  endif()

  if(gridkit_add_executable_TEST)
    add_test(NAME ${target} COMMAND $<TARGET_FILE:${target}>)
  endif()

  install(TARGETS ${target} RUNTIME DESTINATION bin)
endmacro()
