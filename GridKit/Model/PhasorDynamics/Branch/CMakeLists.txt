# [[
#  Author(s):
#    - Cameron Rutherford <cameron.rutherford@pnnl.gov>
#    - Slaven Peles <peless@ornl.gov>
# ]]

set(_install_headers Branch.hpp BranchData.hpp)

if(GRIDKIT_ENABLE_ENZYME)
  gridkit_add_library(
    phasor_dynamics_branch
    SOURCES BranchEnzyme.cpp
    HEADERS ${_install_headers}
    INCLUDE_DIRECTORIES PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include
    LINK_LIBRARIES
      PRIVATE
      ClangEnzymeFlags
      PUBLIC
      GridKit::phasor_dynamics_core
      GridKit::phasor_dynamics_signal
    COMPILE_OPTIONS
      PRIVATE
      -mllvm
      -enzyme-auto-sparsity=1
      -fno-math-errno)
else()
  gridkit_add_library(
    phasor_dynamics_branch
    SOURCES Branch.cpp
    HEADERS ${_install_headers}
    LINK_LIBRARIES PUBLIC GridKit::phasor_dynamics_core GridKit::phasor_dynamics_signal
    INCLUDE_DIRECTORIES PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include)
endif()

gridkit_add_library(
  phasor_dynamics_branch_dependency_tracking
  SOURCES BranchDependencyTracking.cpp
  LINK_LIBRARIES PUBLIC GridKit::phasor_dynamics_core
                 GridKit::phasor_dynamics_signal_dependency_tracking
  INCLUDE_DIRECTORIES PRIVATE ${GRIDKIT_THIRD_PARTY_DIR}/magic-enum/include)

# Link to interface target for all components
target_link_libraries(
  phasor_dynamics_components
  INTERFACE GridKit::phasor_dynamics_branch)
target_link_libraries(
  phasor_dynamics_components_dependency_tracking
  INTERFACE GridKit::phasor_dynamics_branch_dependency_tracking)
