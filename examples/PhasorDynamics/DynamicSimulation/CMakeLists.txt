add_subdirectory(Toy)

if(GRIDKIT_ENABLE_ENZYME)
  add_subdirectory(ACTIVSg10k)
  add_subdirectory(ACTIVSg25kPartitioned)
  add_subdirectory(ACTIVSg70kPartitioned)
endif()
