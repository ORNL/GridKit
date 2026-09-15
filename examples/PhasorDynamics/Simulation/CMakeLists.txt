add_subdirectory(TwoBus)
add_subdirectory(ThreeBus)

if(GRIDKIT_ENABLE_ENZYME)
  add_subdirectory(IEEE39)
  add_subdirectory(ACTIVSg10k)
endif()
