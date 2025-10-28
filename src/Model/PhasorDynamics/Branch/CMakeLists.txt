
# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#    - Slaven Peles <peless@ornl.gov>
# ]]

if(GRIDKIT_ENABLE_ENZYME)
  gridkit_add_library(phasor_dynamics_branch
    SOURCES
      BranchEnzyme.cpp
    LINK_LIBRARIES
      ClangEnzymeFlags
    COMPILE_OPTIONS
      -mllvm -enzyme-auto-sparsity=1 -fno-math-errno)
else()
  gridkit_add_library(phasor_dynamics_branch SOURCES Branch.cpp)
endif()

gridkit_add_library(phasor_dynamics_branch_dependency_tracking
  SOURCES BranchDependencyTracking.cpp)

# Link to interface target for all components
target_link_libraries(gridkit_phasor_dynamics_components
  INTERFACE GRIDKIT::phasor_dynamics_branch)
target_link_libraries(gridkit_phasor_dynamics_components_dependency_tracking
  INTERFACE GRIDKIT::phasor_dynamics_branch_dependency_tracking)
