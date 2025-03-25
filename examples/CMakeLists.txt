#

add_subdirectory(MatPowerTesting)
add_subdirectory(LinearAlgebra)
add_subdirectory(DistributedGeneratorTest)

if(TARGET SUNDIALS::kinsol)
  add_subdirectory(Grid3Bus)
endif()

if(TARGET SUNDIALS::idas)
  add_subdirectory(AdjointSensitivity)
  add_subdirectory(RLCircuit)
  add_subdirectory(Microgrid)
  add_subdirectory(ScaleMicrogrid)
  add_subdirectory(PhasorDynamics)
  if(GRIDKIT_ENABLE_IPOPT)
    add_subdirectory(DynamicConstrainedOpt)
    add_subdirectory(GenConstLoad)
    add_subdirectory(GenInfiniteBus)
    add_subdirectory(ParameterEstimation)
  endif()
endif()

if(GRIDKIT_ENABLE_ENZYME)
  add_subdirectory(Enzyme)
endif()
