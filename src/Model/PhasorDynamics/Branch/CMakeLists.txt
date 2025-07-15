
# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#    - Slaven Peles <peless@ornl.gov>
# ]]

gridkit_add_library(phasor_dynamics_branch
  SOURCES
    Branch.cpp
  OUTPUT_NAME
    gridkit_phasor_dynamics_branch)

# # For below to works CMake min version 3.13 is required
# target_link_libraries(gridkit_phasor_dynamic_components INTERFACE GRIDKIT::phasor_dynamics_branch)
