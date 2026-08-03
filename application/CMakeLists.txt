if(TARGET SUNDIALS::idas)
  add_subdirectory(PhasorDynamics)
endif()

if(GRIDKIT_ENABLE_EIGEN AND TARGET SUNDIALS::idas)
  add_subdirectory(EMT)
endif()
