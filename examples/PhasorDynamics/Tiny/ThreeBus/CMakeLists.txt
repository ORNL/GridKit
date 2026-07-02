add_subdirectory(Basic)
add_subdirectory(Classical)
add_subdirectory(ConstantSource)
add_subdirectory(ZipLoad)

if(GRIDKIT_ENABLE_ARROW)
  add_subdirectory(Arrow)
endif()
