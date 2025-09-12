# [[
#  Author(s):
#    - Slaven Peles <peless@ornl.gov>
#]]

add_subdirectory(DistributedGeneratorTest)

if(TARGET SUNDIALS::ida)
  if(GRIDKIT_ENABLE_SUNDIALS_SPARSE)
    add_subdirectory(RLCircuit)
    add_subdirectory(Microgrid)
    add_subdirectory(ScaleMicrogrid)
  endif()
endif()
