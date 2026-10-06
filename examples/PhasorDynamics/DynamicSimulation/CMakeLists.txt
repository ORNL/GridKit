add_subdirectory(Toy)

if(GRIDKIT_ENABLE_ENZYME)
  add_subdirectory(ACTIVSg10k)
  add_subdirectory(WECC240Relay)
endif()
