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
endif()

if(GRIDKIT_ENABLE_IPOPT AND GRIDKIT_ENABLE_EIGEN)
  add_subdirectory(Optimization)
endif()
