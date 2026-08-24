add_subdirectory(Tiny)
add_subdirectory(Small)

if(GRIDKIT_ENABLE_ENZYME)
  add_subdirectory(Medium)
  add_subdirectory(Large)
  add_subdirectory(Validation)
endif()
