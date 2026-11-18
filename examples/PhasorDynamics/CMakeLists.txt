add_subdirectory(DynamicSimulation)
add_subdirectory(ContingencyAnalysis)
add_subdirectory(Validation)

if(GRIDKIT_ENABLE_ZMQ AND GRIDKIT_ENABLE_ENZYME)
  add_subdirectory(CoSim)
endif()
