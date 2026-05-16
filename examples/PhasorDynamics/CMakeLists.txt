add_subdirectory(Tiny)
add_subdirectory(Small)
add_subdirectory(Validation)

if(GRIDKIT_ENABLE_ENZYME)
  add_subdirectory(Medium)
  add_subdirectory(Large)
endif()
