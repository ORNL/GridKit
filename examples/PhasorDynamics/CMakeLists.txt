add_subdirectory(DynamicSimulation)
add_subdirectory(ContingencyAnalysis)
add_subdirectory(Validation)

if(GRIDKIT_ENABLE_ENZYME AND GRIDKIT_ENABLE_IPOPT)
  add_subdirectory(OptimalDispatch)
endif()
