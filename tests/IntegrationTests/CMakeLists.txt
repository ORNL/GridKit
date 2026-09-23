add_subdirectory(PhasorDynamics)

if(TARGET SUNDIALS::idas
   AND GRIDKIT_ENABLE_ENZYME
   AND GRIDKIT_ENABLE_IPOPT)
  add_subdirectory(OptimalPowerFlow)
endif()
