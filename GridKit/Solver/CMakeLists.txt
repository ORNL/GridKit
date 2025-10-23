# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#    - Slaven Peles <slaven.peles@pnnl.gov>
#]]

if(TARGET SUNDIALS::kinsol)
  add_subdirectory(SteadyState)
endif()

if(TARGET SUNDIALS::idas)
  add_subdirectory(Dynamic)
  if(GRIDKIT_ENABLE_IPOPT)
    add_subdirectory(Optimization)
  endif()
endif()
