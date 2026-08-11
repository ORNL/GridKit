# [[
#  Author(s):
#    - Slaven Peles <peless@ornl.gov>
#]]

add_subdirectory(DistributedGeneratorTest)

if(TARGET SUNDIALS::idas)
  if(GRIDKIT_ENABLE_SUNDIALS_SPARSE)
    add_subdirectory(RLCircuit)
    add_subdirectory(Microgrid)
    add_subdirectory(ScaleMicrogrid)
    add_subdirectory(Partition)
  endif()
endif()
