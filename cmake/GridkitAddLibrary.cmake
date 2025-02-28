
# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#]]

# add_library macro loosely based on https://github.com/LLNL/sundials/blob/master/cmake/macros/SundialsAddLibrary.cmake

macro(gridkit_add_library target)

  set(options STATIC_ONLY SHARED_ONLY)
  set(oneValueArgs OUTPUT_NAME)
  set(multiValueArgs SOURCES LINK_LIBRARIES INCLUDE_DIRECTORIES)

  # parse arguments
  cmake_parse_arguments(gridkit_add_library
    "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  # library types to create
  set(_libtypes "")
  if(GRIDKIT_BUILD_STATIC AND (NOT gridkit_add_library_SHARED_ONLY))
    set(_libtypes "STATIC")
  endif()
  if(GRIDKIT_BUILD_SHARED AND (NOT gridkit_add_library_STATIC_ONLY))
    set(_libtypes "${_libtypes};SHARED")
  endif()

  # build libraries
  foreach(_libtype ${_libtypes})
    # library suffix
    if(${_libtype} MATCHES "STATIC")
      set(_lib_suffix "_static")
    else()
      set(_lib_suffix "_shared")
    endif()

    # source files
    set(sources ${gridkit_add_library_SOURCES})

    # obj target needs a unique name
    set(obj_target ${target}_obj${_lib_suffix})

    # -- Create object library --

    add_library(${obj_target} OBJECT ${sources})
    
    if(gridkit_add_library_LINK_LIBRARIES)
      if(${_lib_type} MATCHES "STATIC")
        append_static_suffix(gridkit_add_library_LINK_LIBRARIES _all_libs)
      else()
        set(_all_libs ${gridkit_add_library_LINK_LIBRARIES})
      endif()
      target_link_libraries(${obj_target} ${_all_libs})
    endif()
    
    # object files going into shared libs need PIC code
    set_target_properties(${obj_target} PROPERTIES POSITION_INDEPENDENT_CODE TRUE)

    # set target name
    set(_actual_target_name ${target}${_lib_suffix})

    add_library(${_actual_target_name} ${_libtype} $<TARGET_OBJECTS:${obj_target}>)

    if(gridkit_add_library_LINK_LIBRARIES)
      target_link_libraries(${_actual_target_name} ${gridkit_add_library_LINK_LIBRARIES})
    endif()

    add_library(GRIDKIT::${target} ALIAS ${_actual_target_name})

    # Set output name
    if(gridkit_add_library_OUTPUT_NAME)
      set_target_properties(${_actual_target_name} PROPERTIES
        OUTPUT_NAME ${gridkit_add_library_OUTPUT_NAME}
        CLEAN_DIRECT_OUTPUT 1)
    else()
      set_target_properties(${_actual_target_name} PROPERTIES
        OUTPUT_NAME ${target}
        CLEAN_DIRECT_OUTPUT 1)
    endif()

    # Set the library version
    set_target_properties(${_actual_target_name} PROPERTIES
      VERSION ${PACKAGE_VERSION}
      SOVERSION ${PACKAGE_VERSION_MAJOR})

    install(TARGETS ${_actual_target_name} DESTINATION lib EXPORT gridkit-targets)
  endforeach()
endmacro()


macro(append_static_suffix libs_in libs_out)
  set(${libs_out} "")
  set(_STATIC_LIB_SUFFIX "_static")
  foreach(_lib ${${libs_in}})
    if(TARGET ${_lib}${_STATIC_LIB_SUFFIX})
      list(APPEND ${libs_out} ${_lib}${_STATIC_LIB_SUFFIX})
    else()
      list(APPEND ${libs_out} ${_lib})
    endif()
  endforeach()
endmacro()
