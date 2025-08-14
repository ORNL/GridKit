
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

gridkit_add_library(phasor_dynamics_branch_dependency_tracking
  SOURCES
  BranchDependencyTracking.cpp
  OUTPUT_NAME
    gridkit_phasor_dynamics_branch_dependency_tracking)

# Link to interface target for all components
target_link_libraries(gridkit_phasor_dynamics_components INTERFACE GRIDKIT::phasor_dynamics_branch)
