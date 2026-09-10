add_subdirectory(TwoBus)
add_subdirectory(ThreeBus)
add_subdirectory(TenGen)

if(GRIDKIT_ENABLE_ENZYME)
  add_subdirectory(ACTIVSg10k)
endif()
