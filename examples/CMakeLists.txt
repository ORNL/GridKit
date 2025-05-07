# [[
#  Author(s):
#    - Slaven Peles <peless@ornl.gov>
#]]

add_subdirectory(MatPowerTesting)
add_subdirectory(LinearAlgebra)

add_subdirectory(PowerFlow)
add_subdirectory(PowerElectronics)
add_subdirectory(Experimental)

if(TARGET SUNDIALS::idas)
  add_subdirectory(PhasorDynamics)
endif()

if(GRIDKIT_ENABLE_ENZYME)
  add_subdirectory(Enzyme)
endif()
