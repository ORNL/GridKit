add_subdirectory(Toy)

if(GRIDKIT_ENABLE_ENZYME)
  add_subdirectory(IEEE39)
  add_subdirectory(ACTIVSg10k)
endif()
