add_subdirectory(DynamicSimulation)
add_subdirectory(ContingencyAnalysis)
add_subdirectory(Validation)

if(GRIDKIT_ENABLE_ZMQ)
  add_subdirectory(CoSimulation)
endif()
