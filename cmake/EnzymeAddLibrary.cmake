# 
#[[

Finds Enzyme Clang plugin

User may set:
- ENZYME_DIR

Author(s):
- Asher Mancinelli <ashermancinelli@gmail.com>
- Nicholson Koukpaizan <koukpaizannk@ornl.gov>

]]

macro(enzyme_build_object)
  set(options)
  set(oneValueArgs NAME)
  set(multiValueArgs SOURCES LINK_LIBRARIES INCLUDE_DIRECTORIES)
  cmake_parse_arguments(enzyme_build_object "${options}" "${oneValueArgs}"
    "${multiValueArgs}" ${ARGN})

  set(PHASE2 "${CMAKE_CURRENT_BINARY_DIR}/${enzyme_build_object_NAME}.bc")
  set(PHASE3 "${CMAKE_CURRENT_BINARY_DIR}/${enzyme_build_object_NAME}_enzyme.ll")
  set(PHASE4 "${CMAKE_CURRENT_BINARY_DIR}/${enzyme_build_object_NAME}_opt.ll")
  set(PHASE5 "${CMAKE_CURRENT_BINARY_DIR}/${enzyme_build_object_NAME}")

  set(OBJS "")
  set(includes "${enzyme_build_object_INCLUDE_DIRECTORIES}")

  foreach(lib ${enzyme_build_object_LINK_LIBRARIES})
    get_target_property(include ${lib} INCLUDE_DIRECTORIES)
    set(includes "${includes}" ${include})

    get_target_property(libsource ${lib} SOURCES)
    string(FIND "${libsource}" "TARGET" found)
    if(NOT(${found} EQUAL -1))
      list(APPEND LINKER_FLAGS "-Wl,${libsource}")
    endif()
  endforeach()

  foreach(dir ${includes})
    if(EXISTS ${dir})
      list(APPEND INCLUDE_COMPILER_LIST "-I${dir}")
    endif()
  endforeach()

  foreach(SRC ${enzyme_build_object_SOURCES})
    set(PHASE0 "${CMAKE_CURRENT_SOURCE_DIR}/${SRC}")
    set(PHASE1 "${CMAKE_CURRENT_BINARY_DIR}/${enzyme_build_object_NAME}_${SRC}_compile.o")
    add_custom_command(
      DEPENDS ${PHASE0} 
      OUTPUT ${PHASE1}
      COMMAND ${CMAKE_CXX_COMPILER} -flto -c ${PHASE0} ${INCLUDE_COMPILER_LIST} -O2 -fno-vectorize -ffast-math -fno-unroll-loops -fpass-plugin=${ENZYME_CLANG_PLUGIN_LIBRARY} -Xclang -load -Xclang ${ENZYME_CLANG_PLUGIN_LIBRARY} -mllvm -enable-load-pre=0 -mllvm -enzyme-auto-sparsity=1 -o ${PHASE1}
      COMMENT "Compiling ${SRC} to object file for target ${enzyme_build_object_NAME}"
      )
    set(OBJS "${OBJS} ${PHASE1}")
  endforeach()

  cmake_language(EVAL CODE "
  add_custom_command(
    DEPENDS ${OBJS}
    OUTPUT ${PHASE2}
    COMMAND ${GRIDKIT_LLVM_LINK} ${OBJS} -o ${PHASE2}
    COMMENT \"Linking object files to LLVM bytecode for target ${enzyme_build_object_NAME}\"
    )
  ")

  add_custom_command(
    DEPENDS ${PHASE2}
    OUTPUT ${PHASE3}
    COMMAND ${GRIDKIT_OPT} ${PHASE2} -load-pass-plugin=${ENZYME_LLVM_PLUGIN_LIBRARY} -passes=enzyme -o ${PHASE3} -S
    COMMENT "Running Enzyme opt pass on target ${enzyme_build_object_NAME}"
    )

  add_custom_command(
    DEPENDS ${PHASE3}
    OUTPUT ${PHASE4}
    COMMAND ${GRIDKIT_OPT} ${PHASE3} -O2 -o ${PHASE4} -S
    COMMENT "Running remaining opt passes on target ${enzyme_build_object_NAME}"
    )

  add_custom_command(
    DEPENDS ${PHASE4}
    OUTPUT ${PHASE5}
    COMMAND ${CMAKE_CXX_COMPILER} -c ${PHASE4} -o ${PHASE5}
    COMMENT "Generating optimized object file for target ${enzyme_build_object_NAME}"
    )
endmacro()

macro(enzyme_add_executable)
  set(options)
  set(oneValueArgs NAME)
  set(multiValueArgs SOURCES LINK_LIBRARIES INCLUDE_DIRECTORIES)
  cmake_parse_arguments(enzyme_add_executable "${options}" "${oneValueArgs}"
    "${multiValueArgs}" ${ARGN})

  enzyme_build_object(
    NAME "${enzyme_add_executable_NAME}.o"
    SOURCES ${enzyme_add_executable_SOURCES}
    LINK_LIBRARIES ${enzyme_add_executable_LINK_LIBRARIES}
    INCLUDE_DIRECTORIES ${enzyme_add_executable_INCLUDE_DIRECTORIES}
  )

  add_executable("${enzyme_add_executable_NAME}" "${enzyme_add_executable_NAME}.o")
  set_target_properties("${enzyme_add_executable_NAME}" PROPERTIES LINKER_LANGUAGE CXX)
  target_link_libraries("${enzyme_add_executable_NAME}" ${enzyme_add_executable_LINK_LIBRARIES})
endmacro()

macro(enzyme_add_library)
  set(options)
  set(oneValueArgs NAME)
  set(multiValueArgs SOURCES LINK_LIBRARIES INCLUDE_DIRECTORIES)
  cmake_parse_arguments(enzyme_add_library "${options}" "${oneValueArgs}"
    "${multiValueArgs}" ${ARGN})

  enzyme_build_object(
    NAME "${enzyme_add_library_NAME}.o"
    SOURCES ${enzyme_add_library_SOURCES}
    LINK_LIBRARIES ${enzyme_add_library_LINK_LIBRARIES}
    INCLUDE_DIRECTORIES ${enzyme_add_library_INCLUDE_DIRECTORIES}
  )

  add_library("${enzyme_add_library_NAME}" "${enzyme_add_library_NAME}.o")
  set_target_properties("${enzyme_add_library_NAME}" PROPERTIES LINKER_LANGUAGE CXX)
  target_link_libraries("${enzyme_add_library_NAME}" ${enzyme_add_library_LINK_LIBRARIES})
endmacro()
