# [[
#  Author(s):
#    - Slaven Peles <peless@ornl.gov>
#]]

add_subdirectory(DistributedGeneratorTest)

if(TARGET SUNDIALS::idas)
  add_subdirectory(RLCircuit)
  add_subdirectory(Microgrid)
  add_subdirectory(ScaleMicrogrid)
endif()
