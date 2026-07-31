add_subdirectory(Toy)

if(GRIDKIT_ENABLE_ENZYME)
  add_subdirectory(ACTIVSg2000)
  add_subdirectory(ACTIVSg10k)
endif()
