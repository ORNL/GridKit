if(TARGET SUNDIALS::idas)
  add_subdirectory(PhasorDynamics)
endif()

if(TARGET model_opf AND TARGET solvers_opt_ipopt)
  add_subdirectory(OPF)
endif()
