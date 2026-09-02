if(TARGET SUNDIALS::idas)
  add_subdirectory(EMT)
  add_subdirectory(PhasorDynamics)
endif()
