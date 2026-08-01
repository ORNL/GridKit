if(TARGET SUNDIALS::idas)
  add_subdirectory(PhasorDynamics)
endif()

if(TARGET solvers_opt)
  add_subdirectory(OPF)
endif()
