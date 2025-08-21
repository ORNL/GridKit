# [[
#  Author(s):
#    - Slaven Peles <peless@ornl.gov>
#]]

if(TARGET GRIDKIT::solvers_dyn)
  add_subdirectory(AdjointSensitivity)
  if(GRIDKIT_ENABLE_IPOPT)
    add_subdirectory(DynamicConstrainedOpt)
    add_subdirectory(GenConstLoad)
    add_subdirectory(GenInfiniteBus)
    add_subdirectory(ParameterEstimation)
  endif()
endif()
