if(TARGET SUNDIALS::idas)
  add_subdirectory(PhasorDynamics)
endif()

if(GRIDKIT_ENABLE_EIGEN AND TARGET SUNDIALS::idas)
  add_subdirectory(EMT)
endif()

if(GRIDKIT_ENABLE_IPOPT AND GRIDKIT_ENABLE_EIGEN)
  add_subdirectory(Fitting)
endif()
