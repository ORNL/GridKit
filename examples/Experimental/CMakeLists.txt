# [[
#  Author(s):
#    - Slaven Peles <peless@ornl.gov>
#]]

if(TARGET SUNDIALS::idas)
  add_subdirectory(AdjointSensitivity)
endif()
