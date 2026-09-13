add_subdirectory(Simulation)

if(GRIDKIT_ENABLE_ENZYME)
  add_subdirectory(Validation)
endif()
