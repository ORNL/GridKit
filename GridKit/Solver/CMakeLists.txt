# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#    - Slaven Peles <slaven.peles@pnnl.gov>
#]]

if(TARGET SUNDIALS::kinsol)
  add_subdirectory(SteadyState)
endif()

add_subdirectory(Dynamic)

if(TARGET SUNDIALS::idas AND GRIDKIT_ENABLE_IPOPT)
  add_subdirectory(Optimization)
endif()
