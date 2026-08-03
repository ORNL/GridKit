add_subdirectory(FrequencyResponse)

if(GRIDKIT_ENABLE_IPOPT)
  add_subdirectory(UniversalLineModel)
endif()
