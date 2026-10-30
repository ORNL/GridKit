add_subdirectory(Toy)

if(GRIDKIT_ENABLE_ENZYME)
  add_subdirectory(PostgreSQL)
  add_subdirectory(ACTIVSg10k)
  add_subdirectory(Texas7kGas)
endif()
