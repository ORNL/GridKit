# [[
#  Author(s):
#    - Slaven Peles <peless@ornl.gov>
#]]

add_subdirectory(MatPowerTesting)

if(TARGET SUNDIALS::kinsol)
  add_subdirectory(Grid3Bus)
endif()
