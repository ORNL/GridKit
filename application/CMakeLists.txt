if(TARGET SUNDIALS::idas)
  add_subdirectory(PhasorDynamics)
endif()

if(GridKit_ENABLE_NOTEBOOKS)
  add_subdirectory(notebooks)
endif()
