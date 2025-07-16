# [[
#  Author(s):
#    - Slaven Peles <peless@ornl.gov>
#]]

add_subdirectory(DistributedGeneratorTest)

if(TARGET SUNDIALS::idas)
  if(TARGET SUNDIALS:KLU)
    add_subdirectory(RLCircuit)
    add_subdirectory(Microgrid)
    add_subdirectory(ScaleMicrogrid)
  endif()
endif()
