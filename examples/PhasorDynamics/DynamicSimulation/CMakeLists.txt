add_subdirectory(Toy)

if(GRIDKIT_ENABLE_ENZYME)
  add_subdirectory(ACTIVSg10k)
  add_subdirectory(ACTIVSg25k)
  add_subdirectory(ACTIVSg70k)
endif()
